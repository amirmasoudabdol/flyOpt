# Simple Makefile to get all nicely compiled and linked

RM=rm -f

# CC=gcc-4.8
# CC=clang
# CFLAGS=-Wall -g -std=c99 -O2
# -fsanitize=address -fno-omit-frame-pointer
# -O1 -fsanitize=address -fno-omit-frame-pointer
# -O2 -fsanitize=address -fno-omit-frame-pointer

SOURCE_FILES=allocate.c evaluate.c init.c local_search.c \
	recombine.c refine.c report.c sort.c ss.c ssTools.c \
	stats.c update.c
OBJS:=$(patsubst %.c,%.o,$(SOURCE_FILES))

SS_INCLUDES=-I../fly -I../utils -I../ss

MAIN=_ss.o

all: $(MAIN)

$(MAIN): $(OBJS)
	# $(CC) $(CFLAGS) $(OBJS) -o $(MAIN) -lm
	# $(CC) $(CCFLAGS) $(SS_INCLUDES)  -c $<

%.o: %.c ss.h rand.h
	$(CC) $(CFLAGS) $(SSFLAGS) $(SS_INCLUDES) -c $<

clean:
	$(RM) $(OBJS)
