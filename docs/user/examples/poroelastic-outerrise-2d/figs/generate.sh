/Users/danieldouglas/software/texlive/bin/universal-darwin/./pdflatex diagrams.tex
/Users/danieldouglas/software/texlive/bin/universal-darwin/./pdfextract 1 1 diagrams.pdf geometry.pdf
/Users/danieldouglas/software/texlive/bin/universal-darwin/./pdfextract 2 2 diagrams.pdf step01-diagram.pdf
/Users/danieldouglas/software/texlive/bin/universal-darwin/./pdfextract 3 3 diagrams.pdf step02-diagram.pdf
/Users/danieldouglas/software/texlive/bin/universal-darwin/./pdfextract 4 4 diagrams.pdf step03-diagram.pdf

for fin in *.pdf; do
  fout1=`echo $fin | sed -e s/\.pdf/1\.svg/g`
  fout=`echo $fin | sed -e s/\.pdf/\.svg/g`
  mutool convert -F svg -o $fout $fin && mv -f $fout1 $fout
  rm -f diagrams?.svg
done 
