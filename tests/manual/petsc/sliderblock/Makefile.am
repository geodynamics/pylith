TESTS = sliderblock

check_PROGRAMS = sliderblock

sliderblock_SOURCES = \
	SliderBlockApp.cc \
	Formulation.cc \
	QuasistaticPrescribedSlip.cc \
	DynamicPrescribedSlip.cc \
	PrescribedSlip.cc \
	DisplacementBC.cc \
	sliderblock.cc

noinst_HEADERS = \
	SliderBlockApp.hh \
	Formulation.hh \
	QuasistaticPrescribedSlip.hh \
	DynamicPrescribedSlip.cc \
	PrescribedSlip.hh \
	DisplacementBC.hh


AM_CPPFLAGS = $(PETSC_CC_INCLUDES)

sliderblock_LDADD = $(PETSC_LIB)


# End of file
