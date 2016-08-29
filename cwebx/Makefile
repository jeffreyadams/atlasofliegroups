# This file is part of the Atlas of Lie Groups ans Representation project

# A minimal Makefile for in-place compilation of a stripped version of cwebx

# The C compiler should be ANSI compatible. If you change this to another
# brand then gcc, make sure you also change the -ansi flag in CFLAGS below to
# one that tells your compiler to expect only strict ANSI/ISO C source code.
# If you don't, then your compiler might squirm about identifiers that it
# predeclares without us having asked it to, and which conflict with the
# identifiers actually used in the CWEBx source, which then fails to compile.
CC = gcc

# The -ansi flag (for gcc) should guard against spurious #include files
# We keep debugging info around to enable the `+d' option of cweave
CFLAGS = -ansi -DDEBUG

# RM is used below in case rm is aliased
RM= /bin/rm -f

##########  You shouldn't have to change anything after this point #######

.SUFFIXES:
.SUFFIXES: .c .o

.c.o:
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

.PHONY: all cweb
all: cweb
cweb: ctanglex cweavex

ctanglex: ctangle.o common.o
	$(CC) $(CFLAGS) -o ctanglex ctangle.o common.o

cweavex: cweave.o common.o
	$(CC) $(CFLAGS) -o cweavex cweave.o common.o

.PHONY: clean install check

clean:
	$(RM) *.o
	$(RM)  ctanglex cweavex

install: ctanglex cweavex
	@echo "The CWEBx executables remain here for the Atlas project"

check:
	@echo 'Did you expect a check from me?'