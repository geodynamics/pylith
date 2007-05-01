## Process this file with automake to produce Makefile.in
# $Id: Makefile.am 2370 2005-09-09 19:49:10Z leif $

SUBDIRS = libpylith3d module utils examples

abs_top_builddir = @abs_top_builddir@
abs_builddir = @abs_builddir@

bin_PROGRAMS = 
if COND_EMBEDDING
    bin_PROGRAMS += pypylith3d mpipypylith3d
    INTERPRETER = $(bindir)/pypylith3d
    noinstINTERPRETER = $(abs_builddir)/pypylith3d
    noinstPYTHONPATH = $(abs_top_builddir)/python
else
    INTERPRETER = $(PYTHON)
    noinstINTERPRETER = $(PYTHON)
    noinstPYTHONPATH = $(abs_top_builddir)/python:$(abs_top_builddir)/pylith3d/module/.libs
endif

# pylith3d
nobase_pyexec_PYTHON = \
	pylith3d/__init__.py \
	pylith3d/KeywordValueParse.py \
	pylith3d/MaterialModel/__init__.py \
	pylith3d/MaterialModel/IsotropicLinearElastic.py \
	pylith3d/MaterialModel/IsotropicLinearGenMaxwellViscoelastic.py \
	pylith3d/MaterialModel/IsotropicLinearMaxwellViscoelastic.py \
	pylith3d/MaterialModel/IsotropicLinearMaxwellViscoelasticESF.py \
	pylith3d/MaterialModel/IsotropicPowerLawMaxwellViscoelastic.py \
	pylith3d/MaterialModel/IsotropicPowerLawMaxwellViscoelasticZT.py \
	pylith3d/MaterialModel/MaterialModel.py \
	pylith3d/Materials.py \
	pylith3d/PyLithApp.py

# pypylith3d (libpylith3d + pylith3dmodule + embedded Python interpreter)
INCLUDES = -I$(top_srcdir)/pylith3d/module $(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)
pypylith3d_SOURCES = pypylith3d.cc
pypylith3d_LDADD = \
	$(top_builddir)/pylith3d/module/libPyLithLibmodule.a \
	$(top_builddir)/pylith3d/module/libPyLithMeshLibmodule.a \
	$(top_builddir)/pylith3d/libpylith3d/libpylith3d.a
pypylith3d$(EXEEXT): $(pypylith3d_OBJECTS) $(pypylith3d_DEPENDENCIES) 
	@rm -f pypylith3d$(EXEEXT)
	$(CXXLINK) $(PYTHON_LDFLAGS) $(PYTHON_LINKFORSHARED) \
		$(pypylith3d_LDFLAGS) $(pypylith3d_OBJECTS) $(pypylith3d_LDADD) \
		$(PYTHON_BLDLIBRARY) \
		$(PYTHON_EGG_LDFLAGS) \
		$(PETSC_FORTRAN_LIB) $(PETSC_LIB) \
		$(FCLIBS) \
		$(PYTHON_LIBS) $(PYTHON_MODLIBS) $(PYTHON_SYSLIBS) \
		$(LIBS) \
		$(PYTHON_LDLAST)
mpipypylith3d_SOURCES = pypylith3d.cc
mpipypylith3d_CXXFLAGS = -DUSE_MPI
mpipypylith3d_LDADD = \
	$(top_builddir)/pylith3d/module/libPyLithLibmodule.a \
	$(top_builddir)/pylith3d/module/libPyLithMeshLibmodule.a \
	$(top_builddir)/pylith3d/libpylith3d/libpylith3d.a
mpipypylith3d$(EXEEXT): $(mpipypylith3d_OBJECTS) $(mpipypylith3d_DEPENDENCIES) 
	@rm -f mpipypylith3d$(EXEEXT)
	$(CXXLINK) $(PYTHON_LDFLAGS) $(PYTHON_LINKFORSHARED) \
		$(mpipypylith3d_LDFLAGS) $(mpipypylith3d_OBJECTS) $(mpipypylith3d_LDADD) \
		$(PYTHON_BLDLIBRARY) \
		$(PYTHON_EGG_LDFLAGS) \
		$(PETSC_FORTRAN_LIB) $(PETSC_LIB) \
		$(FCLIBS) \
		$(PYTHON_LIBS) $(PYTHON_MODLIBS) $(PYTHON_SYSLIBS) \
		$(LIBS) \
		$(PYTHON_LDLAST)

# applications
bin_SCRIPTS = pylith3dapp.py
do_build = sed -e s%[@]INTERPRETER[@]%$(noinstINTERPRETER)%g -e s%[@]PYTHONPATH[@]%$(noinstPYTHONPATH)%g
do_install = sed -e s%[@]INTERPRETER[@]%$(INTERPRETER)%g -e s%[@]PYTHONPATH[@]%$(PYTHONPATH)%g
pylith3dapp.py:  $(srcdir)/pylith3dapp.py.in Makefile
	$(do_build) <  $(srcdir)/pylith3dapp.py.in > $@ || (rm -f $@ && exit 1)
	chmod +x $@
install-binSCRIPTS: $(bin_SCRIPTS)
	@$(NORMAL_INSTALL)
	test -z "$(bindir)" || $(mkdir_p) "$(DESTDIR)$(bindir)"
	@list='$(bin_SCRIPTS)'; for p in $$list; do \
	  if test -f "$$p.in"; then d=; else d="$(srcdir)/"; fi; \
	  if test -f $$d$$p.in; then \
	    f=`echo "$$p" | sed 's|^.*/||;$(transform)'`; \
	    echo " $(do_install) '$$d$$p.in' > '$(DESTDIR)$(bindir)/$$f'"; \
	    $(do_install) "$$d$$p.in" > "$(DESTDIR)$(bindir)/$$f"; \
	    echo " chmod +x '$(DESTDIR)$(bindir)/$$f'"; \
	    chmod +x "$(DESTDIR)$(bindir)/$$f"; \
	  else :; fi; \
	done
EXTRA_DIST = pylith3dapp.py.in
CLEANFILES = $(bin_SCRIPTS)

## end of Makefile.am
