CC	= g++
FLAGS	= -static -ansi -Wall -pedantic -O3 -s -std=c++11
DFLAGS	= -static -ansi -Wall -pedantic -O0 -g -std=c++11
BIN	= ../bin

NAME1	= find_3node_motifs
NAME2	= motif3_network_changes
NAME3	= motif3_network_changes_details

SOURCE1	= $(NAME1).cpp
SOURCE2	= $(NAME2).cpp
SOURCE3 = $(NAME3).cpp

.SUFFIXES: .cpp .o

.cpp.o:; $(CC) -c $(FLAGS) $<

prog:
	mkdir -p $(BIN); $(CC) -o $(BIN)/$(NAME1) $(FLAGS) $(SOURCE1)
	$(CC) -o $(BIN)/$(NAME2) $(FLAGS) $(SOURCE2)
	$(CC) -o $(BIN)/$(NAME3) $(FLAGS) $(SOURCE3)

debug:
	mkdir -p $(BIN); $(CC) -o $(BIN)/debug.$(NAME1) $(DFLAGS) $(SOURCE1)
	$(CC) -o $(BIN)/debug.$(NAME2) $(DFLAGS) $(SOURCE2)
	$(CC) -o $(BIN)/debug.$(NAME3) $(DFLAGS) $(SOURCE3)

clean:
	rm -f $(BIN)/$(NAME1)
	rm -f $(BIN)/debug.$(NAME1)
	rm -f $(BIN)/$(NAME2)
	rm -f $(BIN)/debug.$(NAME2)
	rm -f $(BIN)/$(NAME3)
	rm -f $(BIN)/debug.$(NAME3)
