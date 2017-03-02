s/\\item \[{/\\facilityitem{/g
s/}\] /}{/g
s/^\\facilityitem/\}\n\\facilityitem/g
s/{description}/{inventory}/g
s/~/ /g
s/{lyxcode}/{shell}/g
/^\\begin{centering}/d
/^\\par\\end{centering}/d
/^\\noindent \\begin{center}/d
/^\\par\\end{center}/d
