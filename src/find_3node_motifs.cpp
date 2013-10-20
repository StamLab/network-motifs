/*
  Author: Shane J. Neph
*/


#include <algorithm>
#include <cstddef>
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
  typedef std::map< std::string, std::set<std::string> > NetworkType;
  typedef NetworkType BidirEdges;
  typedef NetworkType UniEdges;

  struct ByLine : public std::string {
    friend std::istream& operator>>(std::istream& is, ByLine& b) {
      std::getline(is, b);
      return(is);
    }
  };

  struct Help {};

  struct Input {
    Input(int argc, char** argv);

    const BidirEdges& AllBidirectionalEdges() const { return allBis_; }
    const BidirEdges& BidirectionalEdges() const { return bis_; }
    const UniEdges& UnidirectionalInputEdges() const { return ins_; }
    const UniEdges& UnidirectionalOutputEdges() const { return outs_; }

    static std::string Usage() {
      std::string msg = "find_3node_motifs <input-graph>";
      msg += "\n\n  <input-graph> is a file with rows that look like:\nA   B\n  where a tab separates the node labels and A->B in your graph";
      return msg;
    }

  private:
    BidirEdges allBis_;
    BidirEdges bis_;
    UniEdges outs_;
    UniEdges ins_;
  };

  void find_motifs(const Input& input);
} // unnamed


//========
// main()
//========
int main(int argc, char** argv) {
  try {
    Input input(argc, argv);
    find_motifs(input);
    return EXIT_SUCCESS;
  } catch(Help& h) {
    std::cout << Input::Usage() << std::endl;
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
  void ffl(const Input& input) {
    const UniEdges& unidirOutEdges = input.UnidirectionalOutputEdges();
    for ( auto& p : unidirOutEdges ) {
      for ( auto& v : p.second ) {
        auto w = unidirOutEdges.find(v);
        if ( w != unidirOutEdges.end() ) {
          // Print node with two outgoing edges in FFL, then edge with 1 in and 1 out, then the sink
          std::vector<std::string> z;
          std::set_intersection(p.second.begin(), p.second.end(),
                                w->second.begin(), w->second.end(),
                                std::back_inserter(z));
          for ( auto& zout : z )
            std::cout << "FFL:\t" << p.first << "\t" << v << "\t" << zout << std::endl;
        }
      } // for
    } // for
  }

  void tre_loop(const Input& input) {
    const UniEdges& unidirOutEdges = input.UnidirectionalOutputEdges();
    const UniEdges& unidirInEdges = input.UnidirectionalInputEdges();
    for ( auto& p : unidirOutEdges ) {
      auto k = unidirInEdges.find(p.first);
      if ( k != unidirInEdges.end() ) {
        for ( auto& v : p.second ) {
          if ( v < p.first ) {
            auto w = unidirOutEdges.find(v);
            if ( w != unidirOutEdges.end() ) {
              std::vector<std::string> z;

              std::set_intersection(w->second.begin(), w->second.end(),
                                    k->second.begin(), k->second.end(),
                                    std::back_inserter(z));
              for ( auto& zout : z ) {
                if ( v < zout )
                  std::cout << "3-Loop:\t" << p.first << "\t" << v << "\t" << zout << std::endl;
              } // for
            }
          }
        } // for
      }
    } // for
  }

  void tre_chain(const Input& input) {
    const BidirEdges& allBidirEdges = input.AllBidirectionalEdges();
    const UniEdges& unidirOutEdges = input.UnidirectionalOutputEdges();
    const UniEdges& unidirInEdges = input.UnidirectionalInputEdges();
    for ( auto& p : unidirOutEdges ) {
      for ( auto& v : p.second ) {
        auto w = unidirOutEdges.find(v);
        if ( w != unidirOutEdges.end() ) {
          std::set<std::string> z(p.second.begin(), p.second.end()), y;
          auto k1 = unidirInEdges.find(p.first);
          if ( k1 != unidirInEdges.end() )
            z.insert(k1->second.begin(), k1->second.end());
          auto k2 = allBidirEdges.find(p.first);
          if ( k2 != allBidirEdges.end() )
            z.insert(k2->second.begin(), k2->second.end());

          std::set_difference(w->second.begin(), w->second.end(),
                              z.begin(), z.end(),
                              std::inserter(y, y.end()));
          for ( auto& yout : y )
            std::cout << "3-Chain:\t" << p.first << "\t" << v << "\t" << yout << std::endl;
        }
      } // for
    } // for
  }

  void v_out(const Input& input) {
    const BidirEdges& allBidirEdges = input.AllBidirectionalEdges();
    const UniEdges& unidirOutEdges = input.UnidirectionalOutputEdges();
    const UniEdges& unidirInEdges = input.UnidirectionalInputEdges();
    for ( auto& p : unidirOutEdges ) {
      std::vector<std::string> z;
      std::set<std::string> s;
      for ( auto& v : p.second ) {
        s.clear(); z.clear();
        auto w = unidirOutEdges.find(v);
        auto x = unidirInEdges.find(v);
        auto y = allBidirEdges.find(v);
        if ( w != unidirOutEdges.end() )
          s.insert(w->second.begin(), w->second.end());
        if ( x != unidirInEdges.end() )
          s.insert(x->second.begin(), x->second.end());
        if ( y != allBidirEdges.end() )
          s.insert(y->second.begin(), y->second.end());

        std::set_difference(p.second.begin(), p.second.end(),
                            s.begin(), s.end(), std::back_inserter(z));
        for ( auto& zout : z ) {
          if ( v < zout )
            std::cout << "V-out:\t" << p.first << "\t" << v << "\t" << zout << std::endl;
        } // for
      } // for
    } // for
  }

  void v_in(const Input& input) {
    const BidirEdges& allBidirEdges = input.AllBidirectionalEdges();
    const UniEdges& unidirOutEdges = input.UnidirectionalOutputEdges();
    const UniEdges& unidirInEdges = input.UnidirectionalInputEdges();
    for ( auto& p : unidirInEdges ) {
      std::vector<std::string> z;
      std::set<std::string> s;
      for ( auto& v : p.second ) {
        s.clear(); z.clear();
        auto w = unidirOutEdges.find(v);
        auto x = unidirInEdges.find(v);
        auto y = allBidirEdges.find(v);
        if ( w != unidirOutEdges.end() )
          s.insert(w->second.begin(), w->second.end());
        if ( x != unidirInEdges.end() )
          s.insert(x->second.begin(), x->second.end());
        if ( y != allBidirEdges.end() )
          s.insert(y->second.begin(), y->second.end());

        std::set_difference(p.second.begin(), p.second.end(),
                            s.begin(), s.end(), std::back_inserter(z));
        for ( auto& zout : z ) {
          if ( v < zout )
            std::cout << "V-in:\t" << p.first << "\t" << v << "\t" << zout << std::endl;
        } // for
      } // for
    } // for
  }

  void regulating_mutual(const Input& input) {
    // print node associated with 2 unidir edges, then remaining edges in alphabetical order
    const BidirEdges& bidirEdges = input.BidirectionalEdges();
    const UniEdges& unidirOutEdges = input.UnidirectionalOutputEdges();
    for ( auto& p : bidirEdges ) {
      auto out_p = unidirOutEdges.find(p.first);
      if ( out_p != unidirOutEdges.end() ) {
        for ( auto& v : p.second ) {
          auto out_v = unidirOutEdges.find(v);
          if ( out_v != unidirOutEdges.end() ) {
            std::vector<std::string> z;
            std::set_intersection(out_p->second.begin(), out_p->second.end(),
                                  out_v->second.begin(), out_v->second.end(),
                                  std::back_inserter(z));
            for ( auto& zout : z )
              std::cout << "Regulating-Mutual:\t" << zout << "\t" << p.first << "\t" << v << std::endl;
          }
        } // for
      }
    } // for
  }

  void regulated_mutual(const Input& input) {
    // Print node associated with 2 unidir edges, then remaining edges in alphabetical order
    const BidirEdges& bidirEdges = input.BidirectionalEdges();
    const UniEdges& unidirInEdges = input.UnidirectionalInputEdges();
    for ( auto& p : bidirEdges ) {
      auto in_p = unidirInEdges.find(p.first);
      if ( in_p != unidirInEdges.end() ) {
        for ( auto& v : p.second ) {
          auto in_v = unidirInEdges.find(v);
          if ( in_v != unidirInEdges.end() ) {
            std::vector<std::string> z;
            std::set_intersection(in_p->second.begin(), in_p->second.end(),
                                  in_v->second.begin(), in_v->second.end(),
                                  std::back_inserter(z));
            for ( auto& zout : z )
              std::cout << "Regulated-Mutual:\t" << zout << "\t" << p.first << "\t" << v << std::endl;
          }
        } // for
      }
    } // for
  }

  void clique(const Input& input) {
    // remember bidirEdges.first has A only if A < B and A<->B (otherwise it has B as .first)
    const BidirEdges& bidirEdges = input.BidirectionalEdges();
    for ( auto& p : bidirEdges ) {
      for ( auto& v : p.second ) {
        auto w = bidirEdges.find(v);
        if ( w != bidirEdges.end() ) {
          std::vector<std::string> z;
          std::set_intersection(p.second.begin(), p.second.end(),
                                w->second.begin(), w->second.end(),
                                std::back_inserter(z));
          for ( auto& zout : z )
            std::cout << "Clique:\t" << p.first << "\t" << v << "\t" << zout << std::endl;
        }
      } // for
    } // for
  }

  void semi_clique(const Input& input) {
    const BidirEdges& allBidirEdges = input.AllBidirectionalEdges();
    const UniEdges& unidirOutEdges = input.UnidirectionalOutputEdges();
    for ( auto& p : allBidirEdges ) {
      for ( auto& v : p.second ) {
        auto w = unidirOutEdges.find(v);
        if ( w != unidirOutEdges.end() ) {
          std::vector<std::string> z;
          std::set_intersection(p.second.begin(), p.second.end(),
                                w->second.begin(), w->second.end(),
                                std::back_inserter(z));
          for ( auto& zout : z )
            std::cout << "Semi-Clique:\t" << p.first << "\t" << v << "\t" << zout << std::endl;
        }
      } // for
    } // for
  }

  void mutual_tre_chain(const Input& input) {
    const BidirEdges& allBidirEdges = input.AllBidirectionalEdges();
    const UniEdges& unidirOutEdges = input.UnidirectionalOutputEdges();
    const UniEdges& unidirInEdges = input.UnidirectionalInputEdges();
    for ( auto& p : allBidirEdges ) {
      auto out_p = unidirOutEdges.find(p.first);
      if ( out_p != unidirOutEdges.end() ) {
        for ( auto& v : p.second ) {
          auto in_v = unidirInEdges.find(v);
          if ( in_v != unidirInEdges.end() ) {
            std::vector<std::string> z;
            std::set_intersection(out_p->second.begin(), out_p->second.end(),
                                  in_v->second.begin(), in_v->second.end(),
                                  std::back_inserter(z));
            for ( auto& zout : z )
              std::cout << "Mutual-And-3-Chain:\t" << p.first << "\t" << zout << "\t" << v << std::endl;
          }
        } // for
      }
    } // for
  }

  void mutual_v(const Input& input) {
    const BidirEdges& allBidirEdges = input.AllBidirectionalEdges();
    const UniEdges& unidirOutEdges = input.UnidirectionalOutputEdges();
    const UniEdges& unidirInEdges = input.UnidirectionalInputEdges();
    for ( auto& p : allBidirEdges ) {
      std::vector<std::string> z;
      std::set<std::string> s;
      for ( auto& v : p.second ) {
        s.clear(); z.clear();
        auto w = unidirOutEdges.find(v);
        auto x = unidirInEdges.find(v);
        auto y = allBidirEdges.find(v);
        if ( w != unidirOutEdges.end() )
          s.insert(w->second.begin(), w->second.end());
        if ( x != unidirInEdges.end() )
          s.insert(x->second.begin(), x->second.end());
        if ( y != allBidirEdges.end() )
          s.insert(y->second.begin(), y->second.end());

        std::set_difference(p.second.begin(), p.second.end(),
                            s.begin(), s.end(), std::back_inserter(z));
        for ( auto& zout : z ) {
          if ( v < zout )
            std::cout << "Mutual-V:\t" << p.first << "\t" << v << "\t" << zout << std::endl;
        } // for
      } // for
    } // for
  }

  void mutual_out(const Input& input) {
    const BidirEdges& allBidirEdges = input.AllBidirectionalEdges();
    const UniEdges& unidirOutEdges = input.UnidirectionalOutputEdges();
    const UniEdges& unidirInEdges = input.UnidirectionalInputEdges();
    for ( auto& p : allBidirEdges ) {
      auto out_p = unidirOutEdges.find(p.first);
      if ( out_p != unidirOutEdges.end() ) {
        std::vector<std::string> z;
        std::set<std::string> s;
        for ( auto& v : p.second ) {
          z.clear(); s.clear();
          auto w = unidirOutEdges.find(v);
          auto x = unidirInEdges.find(v);
          auto y = allBidirEdges.find(v);
          if ( w != unidirOutEdges.end() )
            s.insert(w->second.begin(), w->second.end());
          if ( x != unidirInEdges.end() )
            s.insert(x->second.begin(), x->second.end());
          if ( y != allBidirEdges.end() )
            s.insert(y->second.begin(), y->second.end());

          std::set_difference(out_p->second.begin(), out_p->second.end(),
                              s.begin(), s.end(), std::back_inserter(z));
          for ( auto& zout : z )
            std::cout << "Mutual-Out:\t" << p.first << "\t" << v << "\t" << zout << std::endl;
        } // for
      }
    } // for
  }

  void mutual_in(const Input& input) {
    const BidirEdges& allBidirEdges = input.AllBidirectionalEdges();
    const UniEdges& unidirOutEdges = input.UnidirectionalOutputEdges();
    const UniEdges& unidirInEdges = input.UnidirectionalInputEdges();
    for ( auto& p : allBidirEdges ) {
      auto in_p = unidirInEdges.find(p.first);
      if ( in_p != unidirInEdges.end() ) {
        std::vector<std::string> z;
        std::set<std::string> s;
        for ( auto& v : p.second ) {
          z.clear(); s.clear();
          auto w = unidirOutEdges.find(v);
          auto x = unidirInEdges.find(v);
          auto y = allBidirEdges.find(v);
          if ( w != unidirOutEdges.end() )
            s.insert(w->second.begin(), w->second.end());
          if ( x != unidirInEdges.end() )
            s.insert(x->second.begin(), x->second.end());
          if ( y != allBidirEdges.end() )
            s.insert(y->second.begin(), y->second.end());

          std::set_difference(in_p->second.begin(), in_p->second.end(),
                              s.begin(), s.end(), std::back_inserter(z));
          for ( auto& zout : z )
            std::cout << "Mutual-In:\t" << p.first << "\t" << v << "\t" << zout << std::endl;
        } // for
      }
    } // for
  }

  void find_motifs(const Input& input) {
    ffl(input);
    tre_loop(input);
    tre_chain(input);
    v_out(input);
    v_in(input);
    regulating_mutual(input);
    regulated_mutual(input);
    clique(input);
    semi_clique(input);
    mutual_tre_chain(input);
    mutual_v(input);
    mutual_out(input);
    mutual_in(input);
  }

  Input::Input(int argc, char** argv) {
    for ( int i = 1; i < argc; ++i ) {
      if ( std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h" )
        throw(Help());
    } // for
    if ( argc != 2 )
      throw(Usage());

    std::ifstream infile(argv[1]);
    if ( !infile )
      throw(std::string("Unable to find input file: ") + argv[1]);

    std::string nodeA, nodeB;
    ByLine bline;
    std::string::size_type p;
    std::set<std::string> uniques;
    std::size_t rowid = 0;
    while ( infile >> bline ) {
      ++rowid;
      if ( (p = bline.find('\t')) == std::string::npos ) {
        std::stringstream conv; conv << rowid;
        throw("No tab found at row: " + conv.str());
      }

      if ( !uniques.insert(bline).second ) // not a unique edge
        continue; // otherwise, could add A->B, remove if B->A occurs, and re-add A->B later

      nodeA = bline.substr(0, p); nodeB = bline.substr(++p);
      if ( nodeA == nodeB ) // no self-edges in 3-node motifs
        continue;

      auto iterB = outs_.find(nodeB);
      if ( iterB != outs_.end() && iterB->second.find(nodeA) != iterB->second.end() ) {
        // have A->B and B->A.  Remove edges from UniEdges lists.  Put in BidirEdges list.
        iterB->second.erase(nodeA);
        ins_.find(nodeA)->second.erase(nodeB);

        std::string preferred_source = (nodeA < nodeB) ? nodeA : nodeB;
        std::string preferred_target = (nodeA < nodeB) ? nodeB : nodeA;
        auto iterC = bis_.find(preferred_source);
        auto iterD = allBis_.find(preferred_source), iterE = allBis_.find(preferred_target);
        if ( iterC != bis_.end() )
          iterC->second.insert(preferred_target);            
        else {
          std::set<std::string> s;
          s.insert(preferred_target);
          bis_.insert(std::make_pair(preferred_source, s));
        }

        if ( iterD != allBis_.end() )
          iterD->second.insert(preferred_target);
        else {
          std::set<std::string> s;
          s.insert(preferred_target);
          allBis_.insert(std::make_pair(preferred_source, s));
        }

        if ( iterE != allBis_.end() )
          iterE->second.insert(preferred_source);
        else {
          std::set<std::string> s;
          s.insert(preferred_source);
          allBis_.insert(std::make_pair(preferred_target, s));
        }
      } else {
        auto iterA = outs_.find(nodeA);
        if ( iterA == outs_.end() ) {
          std::set<std::string> s;
          s.insert(nodeB);
          outs_.insert(std::make_pair(nodeA, s));
        } else {
          iterA->second.insert(nodeB);
        }
        auto iterB = ins_.find(nodeB);
        if ( iterB == ins_.end() ) {
          std::set<std::string> s;
          s.insert(nodeA);
          ins_.insert(std::make_pair(nodeB, s));
        } else {
          iterB->second.insert(nodeA);
        }
      }
    } // while
  }
} // unnamed
