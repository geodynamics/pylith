# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#                        (C) 2005  All Rights Reserved
#
#  All worldwide rights reserved.  A license to use, copy, modify and
#  distribute this software for non-commercial research purposes only
#  is hereby granted, provided that this copyright notice and
#  accompanying disclaimer is not modified or removed from the software.
#
#  DISCLAIMER:  The software is distributed "AS IS" without any express
#  or implied warranty, including but not limited to, any implied
#  warranties of merchantability or fitness for a particular purpose
#  or any warranty of non-infringement of any current or pending patent
#  rights.  The authors of the software make no representations about
#  the suitability of this software for any particular purpose.  The
#  entire risk as to the quality and performance of the software is with
#  the user.  Should the software prove defective, the user assumes the
#  cost of all necessary servicing, repair or correction.  In
#  particular, neither Rensselaer Polytechnic Institute, nor the authors
#  of the software are liable for any indirect, special, consequential,
#  or incidental damages related to the software, to the maximum extent
#  the law permits.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

include local.def

PROJECT = lithomop3d
PACKAGE = liblithomop3d

include petsc/default.def

PROJ_SAR = $(BLD_LIBDIR)/$(PACKAGE).$(EXT_SAR)
PROJ_TMPDIR = $(BLD_TMPDIR)/$(PROJECT)/$(PACKAGE)
PROJ_CLEAN += $(PROJ_INCDIR) $(PROJ_SAR)

PROJ_SRCS = \
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
	open_ucd.f \
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
        residu.f \
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
        update_state_cmp.f \
        update_state_drv.f \
        viscos.F \
        winklf.f \
        winklr.F \
	write_bc.f \
	write_connect.f \
	write_cooords.f \
	write_diff.f \
        write_element_info.f \
	write_fuldat.f \
        write_global_info.f \
	write_hist.f \
        write_props.f \
	write_skew.f \
	write_slip.f \
        write_sparse_info.f \
	write_split.f \
	write_split_plot.f \
        write_state_cmp.f \
        write_state_drv.f \
        write_stateout.f \
        write_strscomp.f \
        write_subiter.f \
        write_timdat.f \
        write_ucd_header.f \
        write_ucd_mesh.f \
        write_ucd_node_vals.f \
        write_wink.f \
        write_winkx.f


#--------------------------------------------------------------------------

all: $(PROJ_SAR) export

testit: test

#--------------------------------------------------------------------------
# build the shared object

$(PROJ_SAR): product_dirs $(PROJ_OBJS)
	$(CXX) -o $(PROJ_SAR) $(PROJ_OBJS) $(PROJ_LCXX_FLAGS)

#--------------------------------------------------------------------------
#
export:: export-libraries

test:: show-configuration show-project cpp-configuration fortran-configuration


EXPORT_LIBS = $(PROJ_SAR)


# version
# $Id: Make.mm,v 1.20 2005/04/20 18:59:12 willic3 Exp $

# Generated automatically by MakeMill on Tue Mar  2 17:05:23 2004

# End of file 
