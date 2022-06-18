pdflatex diagrams.tex
pdfextract 1 1 diagrams.pdf step01-diagram.pdf
pdfextract 2 2 diagrams.pdf step02-diagram.pdf
pdfextract 3 3 diagrams.pdf step03-diagram.pdf
pdfextract 4 4 diagrams.pdf step04-diagram.pdf
pdfextract 5 6 diagrams.pdf step06-diagram.pdf

for fin in *.pdf; do
  fout1=`echo $fin | sed -e s/\.pdf/1\.svg/g`
  fout=`echo $fin | sed -e s/\.pdf/\.svg/g`
  mutool convert -F svg -o $fout $fin && mv -f $fout1 $fout
  rm -f diagrams?.svg
done 
