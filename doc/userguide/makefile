#XSLSTYLESHEET_DIR=/usr/share/sgml/docbook/xsl-stylesheets
XSLSTYLESHEET_DIR=/sw/share/xml/xsl/docbook-xsl

pylith_userguide.pdf pdf: pylith_userguide.fo
	fop $< pylith_userguide.pdf

pylith_userguide.ps ps: pylith_userguide.fo
	fop -fo $< -ps pylith_userguide.ps

pylith_userguide.xhtml xhtml:
	xsltproc --xinclude -o pylith_userguide.xhtml ${XSLSTYLESHEET_DIR}/xhtml/docbook.xsl pylith_userguide.xml

pylith_userguide.fo:
	xsltproc --xinclude -o $@ ${XSLSTYLESHEET_DIR}/fo/docbook.xsl pylith_userguide.xml


clean:
	-rm pylith_userguide.fo pylith_userguide.pdf
