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
	mesh2Din3D.txt \
	cube2_ascii.gmv \
	cube2_ascii.pset \
	cube2.txt \
	twohex8.exo \
	twohex8.txt

noinst_TMP = \
	cube2_test.txt \
	mesh2Din3D_test.txt \
	twohex8_test.txt

# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/unittests/pytests/meshio/data
export-data: $(noinst_DATA)
	for f in $(noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(noinst_DATA) \
	$(export_datadir)/$(noinst_TMP)


# End of file 
