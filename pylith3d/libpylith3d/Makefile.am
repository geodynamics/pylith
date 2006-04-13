## Process this file with automake to produce Makefile.in
# $Id: Makefile.am,v 1.1 2005/08/24 22:45:51 leif Exp $

noinst_LIBRARIES =
lib_LTLIBRARIES =

if COND_EMBEDDING
    # static library
    noinst_LIBRARIES += libpylith3d.a
else
    # shared library (libtool)
    lib_LTLIBRARIES += libpylith3d.la
endif

# static library
libpylith3d_a_FFLAGS = $(AM_FFLAGS) # hack for automake
libpylith3d_a_SOURCES = $(sources)

# shared library (libtool)
libpylith3d_la_SOURCES = $(sources)
libpylith3d_la_LDFLAGS =
libpylith3d_la_LIBADD = \
	$(PETSC_LIB)

F77 = $(FC)
AM_FFLAGS = \
	$(FCFLAGS) \
	-I$(top_srcdir)/pylith3d/libpylith3d/includes \
	$(PETSC_INCLUDE) \
	$(MPIINCLUDES)
sources = \
	addfor.f \
	addpr.f \
	addsn.f \
	addstf.F \
	adfldp.f \
	assign_wink.f \
	autoprestr.F \
	bdeld_ss.f \
	bmatrixb.f \
	bmatrixn.f \
	bnmatrix.f \
	bsum.f \
	choldc2.f \
	cholsl.f \
	ckdiag.F \
	cklock.f \
	cmp_stiffsz.f \
	const.f \
	convert_case.f \
	create_id.f \
	cross.f \
	disp.f \
	eforce.f \
	elas_matinit_cmp_ss.F \
	elas_strs_cmp_ss.f \
	elas_strs_mat_cmp_ss.F \
	elastc.F \
	fill.f \
	formdf_ss.f \
	formes_ss.f \
	formf_ss.f \
	formrt.f \
	funcs.f \
	get_initial_stress.f \
	get_units.f \
	getder.f \
	getjac.f \
	getmat.f \
	getshapb.f \
	getshapn.f \
	gload_cmp_ss.f \
	gload_drv.f \
	gravld.f \
	id_split.f \
	idisp.f \
	ifill.f \
	indexx.f \
	infcmp.f \
	infellh.f \
	infelqh.f \
	invar.f \
	iquate.f \
	isort.f \
	iterate.F \
	jacobi.f \
	lcoord.f \
	ldisbc.f \
	ldisp.f \
	ldupdat.f \
	lfit.f \
	lflteq.f \
	lnklst.f \
	load.f \
	loadf.f \
	loadx.f \
	local.f \
	localf.f \
	localx.f \
	makemsr.F \
	mat_1.f \
	mat_2.f \
	mat_3.f \
	mat_4.f \
	mat_5.f \
	mat_6.f \
	mat_7.f \
	mat_8.f \
	mat_9.f \
	mat_10.f \
	mat_11.f \
	mat_12.f \
	mat_13.f \
	mat_14.f \
	mat_15.f \
	mat_16.f \
	mat_17.f \
	mat_18.f \
	mat_19.f \
	mat_20.f \
	matinit_drv.F \
	matmod_def.f \
	meansh.f \
	nchar.f \
	nfind.f \
	nnblnk.f \
	open_append.F \
	open_ucd.F \
	plinhex.f \
	plinpyr.f \
	plintet.f \
	plinwedge.f \
	poldcmp.f \
	pquadhex.f \
	pquadpyr.f \
	pquadtet.f \
	pquadwedge.f \
	preshape.f \
	prestr_matinit_cmp_ss.F \
	presurql.f \
	printd.f \
	printf.f \
	printl.f \
	printv.f \
	prntforc.f \
	pskip.f \
	rdisp.f \
	read_bc.f \
	read_connect.f \
	read_coords.f \
	read_diff.f \
	read_fuldat.f \
	read_hist.f \
	read_skew.f \
	read_slip.f \
	read_split.f \
	read_stateout.f \
	read_timdat.f \
	read_wink.f \
	residu.F \
	rpforc.f \
	rsplit.f \
	rstiff.f \
	rstress.f \
	scan_bc.f \
	scan_connect.f \
	scan_coords.f \
	scan_diff.f \
	scan_fuldat.f \
	scan_hist.f \
	scan_skew.f \
	scan_slip.f \
	scan_split.f \
	scan_timdat.f \
	scan_wink.f \
	scan_winkx.f \
	skclear.f \
	skcomp.f \
	sort_elements.f \
	sort_slip_nodes.f \
	sort_split_nodes.f \
	sprod.f \
	stiff_ss.f \
	stiffld.f \
	stress_drv.f \
	stress_mat_drv.F \
	symmet.f \
	td_matinit_cmp_ss.F \
	td_strs_cmp_ss.f \
	td_strs_mat_cmp_ss.F \
	transp.f \
	try_binio.F \
	unlock.f \
	update_state_cmp.f \
	update_state_drv.f \
	viscos.F \
	winklf.f \
	winklr.F \
	write_bc.F \
	write_connect.F \
	write_coords.F \
	write_diff.F \
	write_element_info.F \
	write_fuldat.F \
	write_global_info.F \
	write_hist.F \
	write_props.F \
	write_skew.F \
	write_slip.F \
	write_sparse_info.F \
	write_split.F \
	write_split_plot.F \
	write_state_cmp.F \
	write_state_drv.F \
	write_stateout.F \
	write_strscomp.F \
	write_subiter.F \
	write_timdat.F \
	write_ucd_header.F \
	write_ucd_mesh.F \
	write_ucd_node_vals.F \
	write_wink.F \
	write_winkx.F \
	includes/binio.h \
	includes/elas_matinit_ext.inc \
	includes/elas_strs_ext.inc \
	includes/elmlbl_def.inc \
	includes/elmlbl_dim.inc \
	includes/getstate_ext.inc \
	includes/gtol_dim.inc \
	includes/implicit.inc \
	includes/labelc_def.inc \
	includes/labelc_dim.inc \
	includes/labeld_def.inc \
	includes/labeld_dim.inc \
	includes/labels_def.inc \
	includes/labels_dim.inc \
	includes/materials.inc \
	includes/mat_prt_ext.inc \
	includes/ncodat_def.inc \
	includes/ncodat_dim.inc \
	includes/nconsts.inc \
	includes/ndimens.inc \
	includes/nforce_def.inc \
	includes/nforce_dim.inc \
	includes/npar_def.inc \
	includes/npar_dim.inc \
	includes/nprint_def.inc \
	includes/nprint_dim.inc \
	includes/nrdimens_dim.inc \
	includes/nshape.inc \
	includes/nsysdat_def.inc \
	includes/nsysdat_dim.inc \
	includes/ntimdat_def.inc \
	includes/ntimdat_dim.inc \
	includes/nunits_def.inc \
	includes/nunits_dim.inc \
	includes/nvisdat_def.inc \
	includes/nvisdat_dim.inc \
	includes/prestr_matinit_ext.inc \
	includes/prscal_dim.inc \
	includes/rconsts.inc \
	includes/rgiter_def.inc \
	includes/rgiter_dim.inc \
	includes/rtimdat_def.inc \
	includes/rtimdat_dim.inc \
	includes/td_matinit_ext.inc \
	includes/td_strs_ext.inc \
	includes/td_strs_mat_ext.inc \
	includes/title_dim.inc \
	includes/update_ext.inc

## end of Makefile.am
