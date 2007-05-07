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

subpackage = feassemble
include $(top_srcdir)/subpackage.am

subpkginclude_HEADERS = \
	ExplicitElasticity.hh \
	ExplicitElasticity.icc \
	ImplicitElasticity.hh \
	ImplicitElasticity.icc \
	Integrator.hh \
	IntegratorElasticity.hh \
	IntegratorExplicit.hh \
	IntegratorImplicit.hh \
	Quadrature.hh \
	Quadrature.icc \
	Quadrature1D.hh \
	Quadrature1D.icc \
	Quadrature1Din2D.hh \
	Quadrature1Din2D.icc \
	Quadrature1Din3D.hh \
	Quadrature1Din3D.icc \
	Quadrature2D.hh \
	Quadrature2D.icc \
	Quadrature2Din3D.hh \
	Quadrature2Din3D.icc \
	Quadrature3D.hh \
	Quadrature3D.icc \
	ParameterManager.hh

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
