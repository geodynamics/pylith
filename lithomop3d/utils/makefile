#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#  Example makefile for TECTON sources
#
#  This makefile was set up for a Macintosh (Mac OSX) using the g77
#  compiler.  If you use a different compiler, change the FCOMPL
#  variable.  You will also probably need to specify a different
#  location for the BLAS and LAPACK libraries (LINALG variable).
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFLAGC= -fno-silent -fno-ugly -fno-f77 -ff90 -fpedantic -Wall -fbackslash -v

CFLAGC= 

# opt   = -g

opt   = -O3 -funroll-loops

load  = -lcc_dynamic

# load  =

FCOMPL= g77

CCOMPL= gcc

INCD= includes

# LINALG=  -L/sw/lib -llapack -lf77blas -lcblas -latlas

LINALG=  -Wl,-framework -Wl,vecLib

.SUFFIXES: .F .f .c .s .o .fil

.F.f:
	$(FCOMPL) -F $<

.f.o:
	$(FCOMPL) -c $(FFLAGC) $(opt) $<

.f.fil:
	$(FCOMPL) -il $(FFLAGC) $<

.s.o:
	as $<

.c.o:
	$(CCOMPL) -c $(CFLAGC) $(opt) $<

lh2qh: lh2qh.o
	$(FCOMPL) $(opt) -o lh2qh lh2qh.o

fixqt: fixqt.o
	$(FCOMPL) $(opt) -o fixqt fixqt.o

readucd: readucd.o
	$(FCOMPL) $(opt) -o readucd readucd.o

tetcmp: tetcmp.o
	$(FCOMPL) $(opt) -o tetcmp tetcmp.o

clean:
	/bin/rm -f *.o
