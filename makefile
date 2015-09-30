DIR_INC = ./include
DIR_SRC = ./src
DIR_OBJ = ./obj
DIR_BIN = ./bin
EXT = cc

SRC = $(wildcard ${DIR_SRC}/*.${EXT})
OBJ = $(patsubst %.${EXT},${DIR_OBJ}/%.o,$(notdir ${SRC})) 

TARGET = hhc

BIN_TARGET = ${DIR_BIN}/${TARGET}

CC = mpicxx
LFLAGS = -lblas -llapack
CFLAGS = -Wall -I${DIR_INC} #-DNDEBUG

${BIN_TARGET}:${OBJ}
	$(CC) $(OBJ)  -o $@ $(LFLAGS)
    
${DIR_OBJ}/%.o:${DIR_SRC}/%.${EXT}
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY:clean
clean:
	rm obj/* bin/*
