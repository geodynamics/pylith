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

TESTS = testpetsc.py

check_SCRIPTS = testpetsc.py

noinst_PYTHON = \
	TestPetscApp.py

noinst_DATA = 

noinst_TMP = 

TESTS_ENVIRONMENT = $(PYTHON)


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests/petsc
export-data: $(noinst_PYTHON) $(noinst_DATA)
	for f in $(noinst_PYTHON) $(noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(noinst_PYTHON) \
	$(export_datadir)/$(noinst_DATA) \
	$(export_datadir)/$(noinst_TMP)


# End of file 
