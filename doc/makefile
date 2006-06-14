XSLSTYLESHEET_DIR=/usr/share/sgml/docbook/xsl-stylesheets-1.69.1-5

pylith_userguide.pdf pdf: pylith_userguide.fo
	fop $< pylith_userguide.pdf

pylith_userguide.html html:
	xsltproc --xinclude -o pylith_userguide.html ${XSLSTYLESHEET_DIR}/html/docbook.xsl pylith_userguide.xml

pylith_userguide.fo:
	xsltproc --xinclude -o $@ ${XSLSTYLESHEET_DIR}/fo/docbook.xsl pylith_userguide.xml


clean:
	-rm pylith_userguide.fo pylith_userguide.pdf
