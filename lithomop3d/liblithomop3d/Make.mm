# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#                        (C) 2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

include local.def

PROJECT = lithomop3d
PACKAGE = liblithomop3d

PROJ_SAR = $(BLD_LIBDIR)/$(PACKAGE).$(EXT_SAR)
PROJ_TMPDIR = $(BLD_TMPDIR)/$(PROJECT)/$(PACKAGE)
PROJ_CLEAN += $(PROJ_INCDIR) $(PROJ_SAR)

PROJ_SRCS = \
        addpr.f \
        ckdiag.f \
        cklock.f \
        const.f \
        elastc.f \
        formdf.f \
        formes.f \
        formf.f \
        formk.f \
        formmat.f \
        gload.f \
        iterate.f \
        ldupdat.f \
        load.f \
        loadf.f \
        loadx.f \
        mathist.f \
        printd.f \
        printf.f \
        printl.f \
        prints.f \
        printv.f \
        prntforc.f \
        residu.f \
        stresn.f \
        unlock.f \
        viscos.f \
        winklf.f \
        winklr.f \
        addstf.f \
        gspre.f \
        isort.f \
        lnklst.f \
        makemsr.f \
        pcginv.f \
	adjid.f \
	convert_case.f \
	get_units.f \
	id_split.f \
	local.f \
	localf.f \
	localx.f \
	nchar.f \
	nfind.f \
	nnblnk.f \
	pskip.f \
	read_bc.f \
	read_connect.f \
	read_coords.f \
	read_diff.f \
	read_fuldat.f \
	read_hist.f \
	read_prestr.f \
	read_prop.f \
	read_skew.f \
	read_slip.f \
	read_split.f \
	read_timdat.f \
	read_traction.f \
	read_wink.f \
	read_winkx.f \
	scan_bc.f \
	scan_connect.f \
	scan_coords.f \
	scan_diff.f \
	scan_fuldat.f \
	scan_hist.f \
	scan_prestr.f \
	scan_prop.f \
	scan_skew.f \
	scan_slip.f \
	scan_split.f \
	scan_timdat.f \
	scan_traction.f \
	scan_wink.f \
	scan_winkx.f \
	write_element_info.f \
	write_global_info.f \
	write_sparse_info.f \
	write_strscomp.f \
	write_subiter.f \
        addfor.f \
        addsn.f \
        adfldp.f \
        bdiff.f \
        disp.f \
        fill.f \
        getmat.f \
        ifill.f \
        indexx.f \
        iquate.f \
        lcoord.f \
        ldisbc.f \
        ldisp.f \
        lflteq.f \
        skclear.f \
        symmet.f \
        transp.f \
        bdeldql.f \
        bmatrixql.f \
        bnmatrxql.f \
        choldc2.f \
        cholsl.f \
        cross.f \
        eforceql.f \
        eqstrsql.f \
        esf.f \
        esfcomp.f \
        formrt.f \
        funcs.f \
        gravldql.f \
        infel.f \
        invar.f \
        jacobi.f \
        lfit.f \
        matflg.f \
        matinit.f \
        matprtb.f \
        meanshql.f \
        poldcmp.f \
        presurql.f \
        rdisp.f \
        rpforc.f \
        rsplit.f \
        rstiff.f \
        rstress.f \
        rtsafe.f \
        shapql.f \
        skcomp.f \
        sprod.f \
        stiffldql.f \
        stiffql.f \
        stresld.f \
        zbrac.f \


#--------------------------------------------------------------------------

all: $(PROJ_SAR) export

#--------------------------------------------------------------------------
# build the shared object

$(PROJ_SAR): product_dirs $(PROJ_OBJS)
	$(CXX) -o $(PROJ_SAR) $(PROJ_OBJS) $(LCXXFLAGS)

#--------------------------------------------------------------------------
#
export:: export-headers export-libraries

EXPORT_HEADERS = \


EXPORT_LIBS = $(PROJ_SAR)


# version
# $Id: Make.mm,v 1.1 2004/04/14 21:18:30 willic3 Exp $

# Generated automatically by MakeMill on Tue Mar  2 17:05:23 2004

# End of file 
