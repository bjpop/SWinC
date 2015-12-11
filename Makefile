CC = gcc
CFLAGS = -Wall
SRC = swinc.c
EXE = swinc
EXE_DEBUG = swinc_debug

$(EXE): $(SRC)
	$(CC) $(CFLAGS) -O2 -o $(EXE) $(SRC)

$(EXE_DEBUG): $(SRC)
	$(CC) $(CFLAGS) -g -o $(EXE_DEBUG) $(SRC)

clean:
	/bin/rm -fr $(EXE) $(EXE_DEBUG)
