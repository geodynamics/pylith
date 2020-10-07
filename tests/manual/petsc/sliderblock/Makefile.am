TESTS = sliderblock

check_PROGRAMS = sliderblock

sliderblock_SOURCES = \
	SliderBlockApp.cc \
	Formulation.cc \
	QuasistaticPrescribedSlip.cc \
	DynamicPrescribedSlip.cc \
	QuasistaticSpontaneousRupture.cc \
	PrescribedSlip.cc \
	DisplacementBC.cc \
	Friction.cc \
	sliderblock.cc

noinst_HEADERS = \
	SliderBlockApp.hh \
	Formulation.hh \
	QuasistaticPrescribedSlip.hh \
	DynamicPrescribedSlip.hh \
	QuasistaticSpontaneousRupture.hh \
	PrescribedSlip.hh \
	Friction.hh \
	DisplacementBC.hh


AM_CPPFLAGS = $(PETSC_CC_INCLUDES)

sliderblock_LDADD = $(PETSC_LIB)


# End of file
