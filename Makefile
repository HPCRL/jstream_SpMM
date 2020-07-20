EXE_DIR = exe
BIN_DIR = bin
SRC_DIR = src
OBJ_DIR = obj
DIRS = $(OBJ_DIR) $(BIN_DIR)

CC = icpc

EXE = $(wildcard $(EXE_DIR)/*.cc)
BIN = $(EXE:$(EXE_DIR)/%.cc=$(BIN_DIR)/%.exe)
SRC = $(wildcard $(SRC_DIR)/*.cc)
OBJ = $(SRC:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)

CPPFLAGS += 
CFLAGS += -Wall -Wno-write-strings -g -std=c++11 -O3 -Ofast  -qopenmp $(EXTRA) -march=native -restrict -mkl
#CFLAGS += -Wall -Wno-write-strings -g -std=c++11 -O0 -qopenmp $(EXTRA) -march=native -restrict -mkl
LDFLAGS += 
LDLIBS += 


.PHONY: all clean

all: | ${DIRS} $(BIN) csb

papi: EXTRA += "-D PAPI"
papi: clean dir 
papi: $(BIN)

${DIRS}:
	mkdir $@

$(BIN_DIR)/%.exe: $(OBJ) $(EXE_DIR)/%.cc
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) $^  -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

csb:
	cd CSB;	make spmm_dall; cd ..


clean: 
	$(RM) -rf $(DIRS)
