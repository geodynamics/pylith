MASTERROOT = pylith_overview

TEXMASTER = $(MASTERROOT).tex

TEXSRCS = 

TEXFILES = $(TEXMASTER) $(TEXSRCS)

AUXFILES = $(TEXMASTER:.tex=.aux)

PDFFILE = $(TEXMASTER:.tex=.pdf)

default: pdf

pdf: $(PDFFILE)

again:
	pdflatex $(TEXMASTER)

$(PDFFILE): $(TEXMASTER)
	pdflatex $(TEXMASTER)

clean:
	rm -f $(AUXFILES) $(PDFFILE) \
		$(MASTERROOT).log $(MASTERROOT).out $(MASTERROOT).nav \
		$(MASTERROOT).snm $(MASTERROOT).toc $(MASTERROOT).vrb *~
	rm -fr auto

