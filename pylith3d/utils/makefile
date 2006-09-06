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
FFLAGC= -Wall -fbackslash -v

CFLAGC= 

opt   = -g

# opt   = -O3 -funroll-loops

# load  = -lcc_dynamic

# load  =

FCOMPL= g95

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

lt2qt: lt2qt.o
	$(FCOMPL) $(opt) -o lt2qt lt2qt.o

fixqt: fixqt.o
	$(FCOMPL) $(opt) -o fixqt fixqt.o

blockrot: blockrot.o
	$(FCOMPL) $(opt) -o blockrot blockrot.o

blockrot2: blockrot2.o
	$(FCOMPL) $(opt) -o blockrot2 blockrot2.o

makeucd: makeucd.o
	$(FCOMPL) $(opt) -o makeucd makeucd.o

pylith2ucd: pylith2ucd.o
	$(FCOMPL) $(opt) -o pylith2ucd pylith2ucd.o

readnetgen: readnetgen.o
	$(FCOMPL) $(opt) -o readnetgen readnetgen.o

readucd: readucd.o
	$(FCOMPL) $(opt) -o readucd readucd.o

readucd2: readucd2.o
	$(FCOMPL) $(opt) -o readucd2 readucd2.o

tetcmp: tetcmp.o
	$(FCOMPL) $(opt) -o tetcmp tetcmp.o

faultcalc: faultcalc.o
	$(FCOMPL) $(opt) -o faultcalc faultcalc.o

faultcalc2: faultcalc2.o
	$(FCOMPL) $(opt) -o faultcalc2 faultcalc2.o

scalecoord: scalecoord.o
	$(FCOMPL) $(opt) -o scalecoord scalecoord.o

scaleucd: scaleucd.o
	$(FCOMPL) $(opt) -o scaleucd scaleucd.o

clean:
	/bin/rm -f *.o
