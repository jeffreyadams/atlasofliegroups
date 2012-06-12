CFLAGS = 
CWEAVE = cweave +d
CTANGLE = ctangle
.SUFFIXES: .tex .dvi .w

.w.tex:
	$(CWEAVE) $*

.tex.dvi:	
	tex $*


.w.dvi:
	$(CWEAVE) $*
	tex $*

.w.c:
	$(CTANGLE) $*

.w.o:
	$(CTANGLE) $*
	$(CC) $(CFLAGS) -c $*.c

.w.dvi:
	make $*.tex
	make $*.dvi

.w.o:
	make $*.c
	make $*.o

treeprint: treeprint.w
	make treeprint.c
	make treeprint	

wc: wc.w
	make wc.c
	make wc

wmerge: wmerge.w
	make wmerge.o
	$(CC) $(CFLAGS) -o wmerge wmerge.o ../common.o


	
