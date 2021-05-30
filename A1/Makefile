###############################################################################
# Makefile for assignment 1, Parallel and Distributed Computing 2020.
###############################################################################

CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm

BIN = stencil

all: $(BIN)

stencil: stencil.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
