/*
  Shane J. Neph
  Univeristy of Washington
  March, 2013
*/

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace {
  enum MotifType {
    Vout = 0, Vin, TreChain, MutualIn, MutualOut, MutualV, FFL, TreLoop, RegulatedMutual,
    RegulatingMutual, MutualAnd3Chain, SemiClique, Clique,

    /* guys with a single node whose position is important */
    Vout_diff, Vin_diff,

    /* guys where every node's position matters */
    TreChain_diff, MutualIn_diff, MutualOut_diff, FFL_diff,

    /* more guys with a single node whose position is important */
    RegulatingMutual_diff, RegulatedMutual_diff, MutualV_diff,

    /* more guys where every node's position matters */
    MutualAnd3Chain_diff, SemiClique_diff,

    /* no 3-node motif in target set */
    NoMatch,

    NumberOfBaseTypes = 13, NumberOfTotalTypes = 25
  };

  typedef std::map< std::string, std::set<std::string> > NetworkType;
  typedef NetworkType BidirEdges;
  typedef NetworkType UniEdges;

  typedef std::string Names; // Nodes C, A, B in a FFL will be represented as ABC for quick lookup
  typedef std::vector<std::string> NodeOrder; // Actual node order; only matters for some motif types
  typedef std::map< Names, std::pair<MotifType, NodeOrder> > NodeLookup;

  typedef std::map< MotifType, std::vector<long> > Counts;


  //========
  // ByLine
  //========
  struct ByLine : public std::string {
    friend std::istream& operator>>(std::istream& is, ByLine& b) {
      std::getline(is, b);
      return(is);
    }
  };

  struct Help {};

  //===========
  // CheckArgs
  //===========
  struct CheckArgs {
    CheckArgs(int argc, char**argv) {
      for ( int i = 1; i < argc; ++i ) {
        if ( std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h" )
          throw(Help());
      } // for
      if ( argc != 3 )
        throw(Usage());

      std::ifstream targetfile(argv[1]);
      if ( !targetfile )
        throw(std::string("Unable to find target file: ") + argv[1]);
      std::ifstream reffile(argv[2]);
      if ( !reffile )
        throw(std::string("Unable to find reference file: ") + argv[2]);

      target_ = argv[1];
      ref_ = argv[2];
    }

    std::string TargetFile() const { return target_; }
    std::string ReferenceFile() const { return ref_; }

    static std::string Usage() {
      std::string msg = "<target-network-file> <reference-network-file>";
      msg += "\nHow do those 3-node circuits found in <reference-network-file> map onto the";
      msg += "\n  same nodes in <target-network-file>?";
      msg += "\n Note that each input files should be the results of running a directed";
      msg += "\n  graph through the find_3node_motifs program.";
      return msg;
    }

  private:
    std::string target_, ref_;
  };


  // Fwd declarations
  void motif_evolution(const NodeLookup& target,
                       const NodeLookup& reference,
                       Counts& counts);
  void read_motifs(const std::string& filename,
                   NodeLookup& lookup);
  void spit_rhymes(Counts& counts);
} // unnamed


//========
// main()
//========
int main(int argc, char** argv) {
  try {
    CheckArgs argcheck(argc, argv);
    NodeLookup target, reference;
    read_motifs(argcheck.TargetFile(), target);
    read_motifs(argcheck.ReferenceFile(), reference);

    Counts counts;
    motif_evolution(target, reference, counts);
    return EXIT_SUCCESS;
  } catch(Help& h) {
    std::cout << CheckArgs::Usage() << std::endl;
    return EXIT_SUCCESS;
  } catch(std::string& s) {
    std::cerr << s << std::endl;
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch(...) {
    std::cerr << "Uknown exception" << std::endl;
  }
  return EXIT_FAILURE;
}


namespace {
  char const* get_name(MotifType m);
  MotifType modified_motif_type(MotifType motifType);

  //===============
  // spit_rhymes()
  //===============
  void spit_rhymes(Counts& counts) {
    // spit out an output matrix

    // header
    std::printf("Motif-Type");
    for ( auto i = Vout; i < NumberOfBaseTypes; i = static_cast<MotifType>(i+1) )
      std::printf("\t%s", get_name(i));
    std::printf("\tNo-Match\tMatched-Variant\n");

    // guts
    MotifType diffType; // for cases when FFLs in both target/ref, but nodes rearranged
    for ( auto i = Vout; i < NumberOfBaseTypes; i = static_cast<MotifType>(i+1) ) {
      std::printf("%s", get_name(i));
      try {
        diffType = modified_motif_type(i);
      } catch(std::string& s) {
        diffType = NoMatch;
      }

      for ( auto j = Vout; j < NumberOfBaseTypes; j = static_cast<MotifType>(j+1) )
        std::printf("\t%ld", counts[i][j]);
      std::printf("\t%ld", counts[i][NoMatch]);

      if ( diffType != NoMatch )
        std::printf("\t%ld", counts[i][diffType]);
      else
        std::printf("\t0");
      std::printf("\n");
    } // for
  }

  //=======================
  // modified_motif_type()
  //=======================
  MotifType modified_motif_type(MotifType motifType) {
    switch (motifType) {
      case FFL:
        return FFL_diff;
      case MutualAnd3Chain:
        return MutualAnd3Chain_diff;
      case MutualIn:
        return MutualIn_diff;
      case TreChain:
        return TreChain_diff;
      case SemiClique:
        return SemiClique_diff;
      case MutualOut:
        return MutualOut_diff;
      case Vout:
        return Vout_diff;
      case Vin:
        return Vin_diff;
      case RegulatingMutual:
        return RegulatingMutual_diff;
      case RegulatedMutual:
        return RegulatedMutual_diff;
      case MutualV:
        return MutualV_diff;
      default:
        throw(std::string("Program Error: modified_motif_type()"));
    };
  }

  //=============
  // get_motif()
  //=============
  MotifType get_motif(const std::string& name) {
    if ( name == "Clique:" )
      return Clique;
    else if ( name == "FFL:")
      return FFL;
    else if ( name == "Mutual-And-3-Chain:" )
      return MutualAnd3Chain;
    else if ( name == "Mutual-In:" )
      return MutualIn;
    else if ( name == "Mutual-Out:" )
      return MutualOut;
    else if ( name == "Mutual-V:" )
      return MutualV;
    else if ( name == "Regulated-Mutual:" )
      return RegulatedMutual;
    else if ( name == "Regulating-Mutual:" )
      return RegulatingMutual;
    else if ( name == "Semi-Clique:" )
      return SemiClique;
    else if ( name == "3-Chain:" )
      return TreChain;
    else if ( name == "3-Loop:" )
      return TreLoop;
    else if ( name == "V-in:" )
      return Vin;
    else if ( name == "V-out:" )
      return Vout;
    else
      throw("Unknown Motif Type: " + name);
  }

  //============
  // get_name()
  //============
  char const* get_name(MotifType m) {
    switch (m) {
      case Clique:
        return "Clique";
      case FFL:
        return "FFL";
      case FFL_diff:
        return "FFL-D";
      case MutualAnd3Chain:
        return "Mutual-And-3-Chain";
      case MutualAnd3Chain_diff:
        return "Mutual-And-3-Chain-D";
      case MutualIn:
        return "Mutual-In";
      case MutualIn_diff:
        return "Mutual-In-D";
      case MutualOut:
        return "Mutual-Out";
      case MutualOut_diff:
        return "Mutual-Out-D";
      case MutualV:
        return "Mutual-V";
      case MutualV_diff:
        return "Mutual-V-D";
      case NoMatch:
        return "No-3-Node-Motif";
      case RegulatedMutual:
        return "Regulated-Mutual";
      case RegulatedMutual_diff:
        return "Regulated-Mutual-D";
      case RegulatingMutual:
        return "Regulating-Mutual";
      case RegulatingMutual_diff:
        return "Regulating-Mutual-D";
      case SemiClique:
        return "Semi-Clique";
      case SemiClique_diff:
        return "Semi-Clique-D";
      case TreChain:
        return "3-Chain";
      case TreChain_diff:
        return "3-Chain-D";
      case TreLoop:
        return "3-Loop";
      case Vin:
        return "V-in";
      case Vin_diff:
        return "V-in-D";
      case Vout:
        return "V-out";
      case Vout_diff:
        return "V-out-D";
      default:
        throw(std::string("Program error: get_name()"));
    };
  }

  //===================
  // motif_evolution()
  //===================
  void motif_evolution(const NodeLookup& target,
                       const NodeLookup& reference,
                       Counts& counts) {

    static const int hardcodeNetworkMotifSize = 3;
    counts.clear();
    for ( auto i = Vout; i < NumberOfBaseTypes; ) {
      counts.insert(std::make_pair(i, std::vector<long>(static_cast<MotifType>(NumberOfTotalTypes), 0)));
      i = static_cast<MotifType>(i+1);
    } // for

    bool same = true;
    for ( auto r : reference ) {
      auto t = target.find(r.first);
      if ( t != target.end() ) {
        if ( t->second.first != r.second.first ) {
          const NodeOrder& ref = r.second.second;
          counts[r.second.first][t->second.first]++;
          std::printf("%s", ref[0].c_str());
          for ( int i = 1; i < hardcodeNetworkMotifSize; ++i )
            std::printf("\t%s", ref[i].c_str());
          std::printf("\t%s\t%s\n", get_name(r.second.first), get_name(t->second.first));
        } else {
          const NodeOrder& ref = r.second.second;
          const NodeOrder& trg = t->second.second;
          MotifType motifType = r.second.first;
          std::printf("%s", ref[0].c_str());
          for ( int i = 1; i < hardcodeNetworkMotifSize; ++i )
            std::printf("\t%s", ref[i].c_str());

          switch (motifType) {
            case FFL: case MutualAnd3Chain: /* all order matters */
            case MutualIn: case MutualOut:
            case TreChain: case SemiClique:
              same = true;
              for ( int i = 0; i < hardcodeNetworkMotifSize; ++i ) {
                if ( ref[i] != trg[i] ) {
                  same = false;
                  break;
                }
              } // for
              if ( same ) {
                counts[motifType][motifType]++;
                std::printf("\t%s\t%s\n", get_name(motifType), get_name(motifType));
              } else {
                counts[motifType][modified_motif_type(motifType)]++;
                std::printf("\t%s\tdiff-%s\n", get_name(motifType), get_name(motifType));
              }
              break;

            /*
            the next non-default case statements require detailed
            knowledge of node order given by the find_3node_motifs program
            turns out it's always a single node, and it's the first given
            */
            case Vout: case Vin: case MutualV:
            case RegulatingMutual: case RegulatedMutual:
              if ( ref[0] == trg[0] ) {
                counts[motifType][motifType]++;
                std::printf("\t%s\t%s\n", get_name(motifType), get_name(motifType));
              } else {
                counts[motifType][modified_motif_type(motifType)]++;
                std::printf("\t%s\tdiff-%s\n", get_name(motifType), get_name(motifType));
              }
              break;

            default:
              counts[motifType][motifType]++;
              std::printf("\t%s\t%s\n", get_name(motifType), get_name(motifType));
          };
        }
      } else {
        counts[r.second.first][NoMatch]++;
        const NodeOrder& ref = r.second.second;
        std::printf("%s", ref[0].c_str());
        for ( int i = 1; i < hardcodeNetworkMotifSize; ++i )
          std::printf("\t%s", ref[i].c_str());
        std::printf("\t%s\tNo-Match\n", get_name(r.second.first));
      }
    } // for
  }

  //===============
  // read_motifs()
  //===============
  void read_motifs(const std::string& filename, NodeLookup& lookup) {
    std::ifstream infile(filename.c_str());
    ByLine bline;
    std::size_t linecntr = 0;
    while ( infile >> bline ) {
      std::stringstream lineNum; lineNum << ++linecntr;
      std::string::size_type p, q;
      if ( bline.find("\t\t") != std::string::npos )
        throw("Consecutive tabs found at line: " + lineNum.str() + " in " + filename);
      else if ( bline.find(" ") != std::string::npos )
        throw("Found 1 or more spaces at line: " + lineNum.str() + " in " + filename);

      p = bline.find("\t");
      if ( p == std::string::npos )
        throw("No tabs at line: " + lineNum.str() + " in " + filename);
      MotifType motifType = get_motif(bline.substr(0, p));
      q = bline.find("\t", ++p);
      if ( q == std::string::npos )
        throw("No tab after second column at line: " + lineNum.str() + " in " + filename);
      std::string firstNode = bline.substr(p, q-p);
      p = bline.find("\t", ++q);
      if ( p == std::string::npos )
        throw("No tab after 3rd column at line: " + lineNum.str() + " in " + filename);
      std::string secondNode = bline.substr(q, p-q);
      std::string thirdNode = bline.substr(++p);
      if ( thirdNode.empty() )
        throw("Problem with 3rd column at line: " + lineNum.str() + " in " + filename);
      else if ( thirdNode.find_first_of("\t") != std::string::npos )
        throw("Tab found in what should be the 3rd column at line: " + lineNum.str() + " in " + filename);

      if ( (firstNode == secondNode) || (firstNode == thirdNode) || (secondNode == thirdNode) )
        throw("The same node found more than once in a 3-motif at line: " + lineNum.str() + " in " + filename);

      std::string a = firstNode, b = secondNode, c = thirdNode;
      if ( b > c )
        std::swap(b, c);

      if ( a > b ) {
        std::swap(a, b);
        if ( b > c )
          std::swap(b, c);
      }

      NodeOrder nodes;
      nodes.push_back(firstNode); nodes.push_back(secondNode); nodes.push_back(thirdNode);
      if ( !lookup.insert(std::make_pair(a+b+c, std::make_pair(motifType, nodes))).second )
        throw(std::string("Multiple rows have the same nodes. One is at line: " + lineNum.str() + " in " + filename));
    } // while
  }
} // unnamed
