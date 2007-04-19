# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

noinst_DATA = \
	meshLine.txt \
	meshTri3A.txt \
	meshQuad4A.txt \
	meshTet4A.txt \
	meshHex8A.txt

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/unittests/libtests/faults/data
export-data: $(noinst_DATA)
	for f in $(noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(noinst_DATA) \
	$(export_datadir)/$(noinst_TMP)


# End of file 
