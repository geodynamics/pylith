## Process this file with automake to produce Makefile.in
# $Id: Makefile.am 2370 2005-09-09 19:49:10Z leif $

fast_tests = \
	linhex/slip1/test.sh \
	linhex/timeslip1/test.sh \
	lintet/bm1b/test.sh

slow_tests = \
	linhex/bm1/test.sh \
	linhex/bm1a/test.sh \
	linhex/bm2a/test.sh \
	$(linhex_ps1_tests) \
	lintet/bm1a/test.sh \
	lintet/bm2a/test.sh \
	lintet/bm5/test.sh \
	quadhex/bm1a/test.sh \
	quadtet/bm1a/test.sh

TESTS = \
	$(fast_tests) \
	$(slow_tests)

XFAIL_TESTS = \
	$(linhex_ps1_tests) \
	quadhex/bm1a/test.sh \
	quadtet/bm1a/test.sh

TESTS_ENVIRONMENT = PATH=$(top_builddir)/lithomop3d:$$PATH PYTHONPATH=$(top_srcdir)/lithomop3d:$$PYTHONPATH

CLEANFILES = \
	*/*/*.ascii \
	*/*/*.mesh.* \
	*/*/*.gmesh.* \
	$(linhex_ps1_tests)

$(linhex_ps1_tests): linhex/ps1/test.sh.in
	fileRoot=`echo $@ | sed -e "s|\.sh|/ps1|"`; \
	do_subst="sed -e s|[@]fileRoot[@]|$$fileRoot|g"; \
	$$do_subst < linhex/ps1/test.sh.in > $@ || (rm -f $@ && exit 1)
	chmod +x $@

check: clean check-am
	$(MAKE) $(AM_MAKEFLAGS) clean

quickcheck:
	$(MAKE) $(AM_MAKEFLAGS) TESTS="$(fast_tests)" check

linhex_ps1_tests = \
	linhex/ps1/elas-test1.sh \
	linhex/ps1/elas-test2.sh \
	linhex/ps1/elas-test3.sh \
	linhex/ps1/elas-test4.sh \
	linhex/ps1/elas-test5.sh \
	linhex/ps1/elas-test6.sh \
	linhex/ps1/elas-test7.sh \
	linhex/ps1/elas-test8.sh \
	linhex/ps1/elas-test9.sh \
	linhex/ps1/elas-test10.sh \
	linhex/ps1/elas-test11.sh \
	linhex/ps1/elas-test12.sh \
	linhex/ps1/vis-test1.sh \
	linhex/ps1/vis-test2.sh \
	linhex/ps1/vis-test3.sh \
	linhex/ps1/vis-test4.sh \
	linhex/ps1/vis-test5.sh \
	linhex/ps1/vis-test6.sh \
	linhex/ps1/vis-test7.sh \
	linhex/ps1/vis-test8.sh \
	linhex/ps1/vis-test9.sh \
	linhex/ps1/vis-test10.sh \
	linhex/ps1/vis-test11.sh \
	linhex/ps1/vis-test12.sh

EXTRA_DIST = \
	bm1/bm1.bc \
	bm1/bm1.connect \
	bm1/bm1.coord \
	bm1/bm1.fuldat \
	bm1/bm1.hist \
	bm1/bm1.keyval \
	bm1/bm1.prop \
	bm1/bm1.time \
	bm1/README \
	bm2/bm2.bc \
	bm2/bm2.keyval \
	bm2/bm2.time \
	bm2.com \
	bm2/README \
	linhex/bm1a/bm1a.bc \
	linhex/bm1a/bm1a.connect \
	linhex/bm1a/bm1a.coord \
	linhex/bm1a/bm1a.fuldat \
	linhex/bm1a/bm1a.hist \
	linhex/bm1a/bm1a.keyval \
	linhex/bm1a/bm1a.prop \
	linhex/bm1a/bm1a.statevar \
	linhex/bm1a/bm1a.time \
	linhex/bm1a/test.sh \
	linhex/bm1/bm1.bc \
	linhex/bm1/bm1.connect \
	linhex/bm1/bm1.coord \
	linhex/bm1/bm1.fuldat \
	linhex/bm1/bm1.hist \
	linhex/bm1/bm1.keyval \
	linhex/bm1/bm1.prop \
	linhex/bm1/bm1.statevar \
	linhex/bm1/bm1.time \
	linhex/bm1/test.sh \
	linhex/bm2a/bm2a.bc \
	linhex/bm2a/bm2a.keyval \
	linhex/bm2a/bm2a.time \
	linhex/bm2a/test.sh \
	linhex/patchtest/pt10.bc \
	linhex/patchtest/pt11.bc \
	linhex/patchtest/pt12.bc \
	linhex/patchtest/pt1a.bc \
	linhex/patchtest/pt1.bc \
	linhex/patchtest/pt1.connect \
	linhex/patchtest/pt1.coord \
	linhex/patchtest/pt1.keyval \
	linhex/patchtest/pt1.prop \
	linhex/patchtest/pt1.statevar \
	linhex/patchtest/pt1.time \
	linhex/patchtest/pt2.bc \
	linhex/patchtest/pt3.bc \
	linhex/patchtest/pt4.bc \
	linhex/patchtest/pt5.bc \
	linhex/patchtest/pt6.bc \
	linhex/patchtest/pt7.bc \
	linhex/patchtest/pt8.bc \
	linhex/patchtest/pt9.bc \
	linhex/patchtest/pt-bbar.keyval \
	linhex/patchtest/pt-incomp.prop \
	linhex/patchtest/pt-red.keyval \
	linhex/patchtest/runtests-bbar-comp.com \
	linhex/patchtest/runtests-bbar-incomp.com \
	linhex/patchtest/runtests-full-comp.com \
	linhex/patchtest/runtests-full-incomp.com \
	linhex/patchtest/runtests-red-comp.com \
	linhex/patchtest/runtests-red-incomp.com \
	linhex/ps1/elas-test10/ps1.keyval \
	linhex/ps1/elas-test10/ps1.prop \
	linhex/ps1/elas-test10/ps1.time \
	linhex/ps1/elas-test11/ps1.keyval \
	linhex/ps1/elas-test11/ps1.prop \
	linhex/ps1/elas-test11/ps1.time \
	linhex/ps1/elas-test12/ps1.keyval \
	linhex/ps1/elas-test12/ps1.prop \
	linhex/ps1/elas-test12/ps1.time \
	linhex/ps1/elas-test1/ps1.keyval \
	linhex/ps1/elas-test1/ps1.prop \
	linhex/ps1/elas-test1/ps1.time \
	linhex/ps1/elas-test2/ps1.keyval \
	linhex/ps1/elas-test2/ps1.prop \
	linhex/ps1/elas-test2/ps1.time \
	linhex/ps1/elas-test3/ps1.keyval \
	linhex/ps1/elas-test3/ps1.prop \
	linhex/ps1/elas-test3/ps1.time \
	linhex/ps1/elas-test4/ps1.keyval \
	linhex/ps1/elas-test4/ps1.prop \
	linhex/ps1/elas-test4/ps1.time \
	linhex/ps1/elas-test5/ps1.keyval \
	linhex/ps1/elas-test5/ps1.prop \
	linhex/ps1/elas-test5/ps1.time \
	linhex/ps1/elas-test6/ps1.keyval \
	linhex/ps1/elas-test6/ps1.prop \
	linhex/ps1/elas-test6/ps1.time \
	linhex/ps1/elas-test7/ps1.keyval \
	linhex/ps1/elas-test7/ps1.prop \
	linhex/ps1/elas-test7/ps1.time \
	linhex/ps1/elas-test8/ps1.keyval \
	linhex/ps1/elas-test8/ps1.prop \
	linhex/ps1/elas-test8/ps1.time \
	linhex/ps1/elas-test9/ps1.keyval \
	linhex/ps1/elas-test9/ps1.prop \
	linhex/ps1/elas-test9/ps1.time \
	linhex/ps1/ps1.bc \
	linhex/ps1/ps1.connect \
	linhex/ps1/ps1.coord \
	linhex/ps1/ps1.fuldat \
	linhex/ps1/ps1.statevar \
	linhex/ps1/README \
	linhex/ps1/README.elas-tests \
	linhex/ps1/README.vis-tests \
	linhex/ps1/test.sh.in \
	linhex/ps1/vis-test10/ps1.keyval \
	linhex/ps1/vis-test10/ps1.prop \
	linhex/ps1/vis-test10/ps1.time \
	linhex/ps1/vis-test11/ps1.keyval \
	linhex/ps1/vis-test11/ps1.prop \
	linhex/ps1/vis-test11/ps1.time \
	linhex/ps1/vis-test12/ps1.keyval \
	linhex/ps1/vis-test12/ps1.prop \
	linhex/ps1/vis-test12/ps1.time \
	linhex/ps1/vis-test1/ps1.keyval \
	linhex/ps1/vis-test1/ps1.prop \
	linhex/ps1/vis-test1/ps1.time \
	linhex/ps1/vis-test2/ps1.keyval \
	linhex/ps1/vis-test2/ps1.prop \
	linhex/ps1/vis-test2/ps1.time \
	linhex/ps1/vis-test3/ps1.keyval \
	linhex/ps1/vis-test3/ps1.prop \
	linhex/ps1/vis-test3/ps1.time \
	linhex/ps1/vis-test4/ps1.keyval \
	linhex/ps1/vis-test4/ps1.prop \
	linhex/ps1/vis-test4/ps1.time \
	linhex/ps1/vis-test5/ps1.keyval \
	linhex/ps1/vis-test5/ps1.prop \
	linhex/ps1/vis-test5/ps1.time \
	linhex/ps1/vis-test6/ps1.keyval \
	linhex/ps1/vis-test6/ps1.prop \
	linhex/ps1/vis-test6/ps1.time \
	linhex/ps1/vis-test7/ps1.keyval \
	linhex/ps1/vis-test7/ps1.prop \
	linhex/ps1/vis-test7/ps1.time \
	linhex/ps1/vis-test8/ps1.keyval \
	linhex/ps1/vis-test8/ps1.prop \
	linhex/ps1/vis-test8/ps1.time \
	linhex/ps1/vis-test9/ps1.keyval \
	linhex/ps1/vis-test9/ps1.prop \
	linhex/ps1/vis-test9/ps1.time \
	linhex/slip1/slip1.bc \
	linhex/slip1/slip1.connect \
	linhex/slip1/slip1.coord \
	linhex/slip1/slip1.fuldat \
	linhex/slip1/slip1.keyval \
	linhex/slip1/slip1.prop \
	linhex/slip1/slip1.slip \
	linhex/slip1/slip1.statevar \
	linhex/slip1/slip1.time \
	linhex/slip1/slip1.winkx \
	linhex/slip1/test.sh \
	linhex/timeslip1/test.sh \
	linhex/timeslip1/timeslip1.bc \
	linhex/timeslip1/timeslip1.connect \
	linhex/timeslip1/timeslip1.coord \
	linhex/timeslip1/timeslip1.fuldat \
	linhex/timeslip1/timeslip1.hist \
	linhex/timeslip1/timeslip1.keyval \
	linhex/timeslip1/timeslip1.prop \
	linhex/timeslip1/timeslip1.slip \
	linhex/timeslip1/timeslip1.statevar \
	linhex/timeslip1/timeslip1.time \
	linhex/timeslip1/timeslip1.winkx \
	lintet/bm1a/bm1a.bc \
	lintet/bm1a/bm1a.connect \
	lintet/bm1a/bm1a.coord \
	lintet/bm1a/bm1a.fuldat \
	lintet/bm1a/bm1a.hist \
	lintet/bm1a/bm1a.keyval \
	lintet/bm1a/bm1a.prop \
	lintet/bm1a/bm1a.statevar \
	lintet/bm1a/bm1a.time \
	lintet/bm1a/test.sh \
	lintet/bm1b/bm1b.bc \
	lintet/bm1b/bm1b.connect \
	lintet/bm1b/bm1b.coord \
	lintet/bm1b/bm1b.fuldat \
	lintet/bm1b/bm1b.hist \
	lintet/bm1b/bm1b.keyval \
	lintet/bm1b/bm1b.prop \
	lintet/bm1b/bm1b.statevar \
	lintet/bm1b/bm1b.time \
	lintet/bm1b/test.sh \
	lintet/bm2a/bm2a.bc \
	lintet/bm2a/bm2a.keyval \
	lintet/bm2a/bm2a.time \
	lintet/bm2a/test.sh \
	lintet/bm5/bm5.bc \
	lintet/bm5/bm5.connect \
	lintet/bm5/bm5.coord \
	lintet/bm5/bm5.fuldat \
	lintet/bm5/bm5.keyval \
	lintet/bm5/bm5.prop \
	lintet/bm5/bm5.split \
	lintet/bm5/bm5.statevar \
	lintet/bm5/bm5.time \
	lintet/bm5/test.sh \
	lintet/patchtest/pt10.bc \
	lintet/patchtest/pt11.bc \
	lintet/patchtest/pt12.bc \
	lintet/patchtest/pt1.bc \
	lintet/patchtest/pt1.connect \
	lintet/patchtest/pt1.coord \
	lintet/patchtest/pt1.keyval \
	lintet/patchtest/pt1.prop \
	lintet/patchtest/pt1.statevar \
	lintet/patchtest/pt1.time \
	lintet/patchtest/pt2.bc \
	lintet/patchtest/pt3.bc \
	lintet/patchtest/pt4.bc \
	lintet/patchtest/pt5.bc \
	lintet/patchtest/pt6.bc \
	lintet/patchtest/pt7.bc \
	lintet/patchtest/pt8.bc \
	lintet/patchtest/pt9.bc \
	lintet/patchtest/pt-bbar.keyval \
	lintet/patchtest/pt-incomp.prop \
	lintet/patchtest/pt-red.keyval \
	lintet/patchtest/runtests-bbar-comp.com \
	lintet/patchtest/runtests-bbar-incomp.com \
	lintet/patchtest/runtests-full-comp.com \
	lintet/patchtest/runtests-full-incomp.com \
	lintet/patchtest/runtests-red-comp.com \
	lintet/patchtest/runtests-red-incomp.com \
	quadhex/bm1a/bm1a.bc \
	quadhex/bm1a/bm1a.connect \
	quadhex/bm1a/bm1a.coord \
	quadhex/bm1a/bm1a.fuldat \
	quadhex/bm1a/bm1a.hist \
	quadhex/bm1a/bm1a.keyval \
	quadhex/bm1a/bm1a.prop \
	quadhex/bm1a/bm1a.statevar \
	quadhex/bm1a/bm1a.time \
	quadhex/bm1a/test.sh \
	quadhex/patchtest/pt10.bc \
	quadhex/patchtest/pt11.bc \
	quadhex/patchtest/pt12.bc \
	quadhex/patchtest/pt1.bc \
	quadhex/patchtest/pt1.connect \
	quadhex/patchtest/pt1.coord \
	quadhex/patchtest/pt1.keyval \
	quadhex/patchtest/pt1.prop \
	quadhex/patchtest/pt1.statevar \
	quadhex/patchtest/pt1.time \
	quadhex/patchtest/pt2.bc \
	quadhex/patchtest/pt3.bc \
	quadhex/patchtest/pt4.bc \
	quadhex/patchtest/pt5.bc \
	quadhex/patchtest/pt6.bc \
	quadhex/patchtest/pt7.bc \
	quadhex/patchtest/pt8.bc \
	quadhex/patchtest/pt9.bc \
	quadhex/patchtest/pt-bbar.keyval \
	quadhex/patchtest/pt-incomp.prop \
	quadhex/patchtest/pt-red.keyval \
	quadhex/patchtest/runtests-bbar-comp.com \
	quadhex/patchtest/runtests-bbar-incomp.com \
	quadhex/patchtest/runtests-full-comp.com \
	quadhex/patchtest/runtests-full-incomp.com \
	quadhex/patchtest/runtests-red-comp.com \
	quadhex/patchtest/runtests-red-incomp.com \
	quadtet/bm1a/bm1a.bc \
	quadtet/bm1a/bm1a.connect \
	quadtet/bm1a/bm1a.coord \
	quadtet/bm1a/bm1a.fuldat \
	quadtet/bm1a/bm1a.hist \
	quadtet/bm1a/bm1a.keyval \
	quadtet/bm1a/bm1a-old.bc \
	quadtet/bm1a/bm1a-old.connect \
	quadtet/bm1a/bm1a-old.coord \
	quadtet/bm1a/bm1a.prop \
	quadtet/bm1a/bm1a.statevar \
	quadtet/bm1a/bm1a.time \
	quadtet/bm1a/test.sh \
	quadtet/patchtest/pt10.bc \
	quadtet/patchtest/pt11.bc \
	quadtet/patchtest/pt12.bc \
	quadtet/patchtest/pt1.bc \
	quadtet/patchtest/pt1.connect \
	quadtet/patchtest/pt1.coord \
	quadtet/patchtest/pt1.keyval \
	quadtet/patchtest/pt1.prop \
	quadtet/patchtest/pt1.statevar \
	quadtet/patchtest/pt1.time \
	quadtet/patchtest/pt2.bc \
	quadtet/patchtest/pt3.bc \
	quadtet/patchtest/pt4.bc \
	quadtet/patchtest/pt5.bc \
	quadtet/patchtest/pt6.bc \
	quadtet/patchtest/pt7.bc \
	quadtet/patchtest/pt8.bc \
	quadtet/patchtest/pt9.bc \
	quadtet/patchtest/pt-bbar.keyval \
	quadtet/patchtest/pt-incomp.prop \
	quadtet/patchtest/pt-red.keyval \
	quadtet/patchtest/runtests-bbar-comp.com \
	quadtet/patchtest/runtests-bbar-incomp.com \
	quadtet/patchtest/runtests-full-comp.com \
	quadtet/patchtest/runtests-full-incomp.com \
	quadtet/patchtest/runtests-red-comp.com \
	quadtet/patchtest/runtests-red-incomp.com

## end of Makefile.am
