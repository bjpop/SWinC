CC = gcc
CFLAGS = -Wall
SRC = swinc.c
EXE = swinc
EXE_DEBUG = swinc_debug

main: $(SRC)
	$(CC) $(CFLAGS) -O2 -o $(EXE) $(SRC)


debug: $(SRC)
	$(CC) $(CFLAGS) -g -o $(EXE_DEBUG) $(SRC)


clean:
	/bin/rm -fr $(EXE) $(EXE_DEBUG)
