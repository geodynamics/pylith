#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
#
#  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
#  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
#  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


cdef extern from "config.h":
    pass  # FC_FUNC and FC_FUNC_

cdef extern from "stddef.h":
    ctypedef unsigned int size_t

cimport petsc


cdef extern void assign_wink "FC_FUNC_(assign_wink, ASSIGN_WINK)" (
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *) except *

cdef extern void autoprestr "FC_FUNC(autoprestr, AUTOPRESTR)" (
    petsc.Mat *,      # sparse
    petsc.Vec *,
    petsc.Vec *,
    double *,          # force
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    int *,
    double *,
    double *,          # global
    double *,
    double *,
    double *,
    int *,
    int *,
    double *,
    int *,
    int *,
    int *,            # bc
    double *,
    double *,          # slip
    double *,
    double *,
    double *,
    int *,
    int *,
    double *,
    int *,
    int *,
    int *,
    double *,          # split
    int *,
    double *,
    double *,
    double *,          # stiff
    double *,
    double *,          # element
    double *,
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,            # traction
    double *,
    double *,
    double *,
    int *,
    double *,          # material
    int *,
    double *,          # element type
    double *,
    double *,
    int *,
    double *,          # time data
    double *,
    int *,
    int *,
    int *,
    double *,
    double *,
    int *,
    int *,
    int *,
    double *,
    double *,
    double *,
    int *,
    double *,          # iterations
    double *,          # skew
    int *,            # i/o info
    int *,
    int *,
    int *,
    int *,
    char *,           # files
    char *,
    char *,
    int *,            # PETSc logging
    int *,
    int *,            # error codes
    char *,
    size_t,           # string lengths
    size_t,
    size_t,
    size_t) except *

cdef extern void cmp_stiffsz "FC_FUNC_(cmp_stiffsz, CMP_STIFFSZ)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void create_id "FC_FUNC_(create_id, CREATE_ID)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *) except *

cdef extern void elastc "FC_FUNC(elastc, ELASTC)" (
    petsc.Mat *,       # sparse
    petsc.Vec *,
    petsc.Vec *,
    double *,           # force
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    int *,
    double *,
    double *,           # global
    double *,
    double *,
    double *,
    int *,
    int *,
    double *,
    int *,
    int *,
    int *,             # bc
    double *,
    double *,           # slip
    double *,
    double *,
    double *,
    int *,
    int *,
    double *,
    int *,
    int *,
    int *,
    double *,           # split
    int *,
    double *,
    double *,
    double *,           # stiff
    double *,
    double *,           # element
    double *,
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,             # traction
    double *,
    double *,
    double *,
    int *,
    double *,           # material
    int *,
    double *,           # element type
    double *,
    double *,
    int *,
    double *,           # time data
    double *,
    int *,
    int *,
    int *,
    double *,
    double *,
    int *,
    int *,
    int *,
    double *,
    double *,
    double *,
    int *,
    double *,           # iterations
    double *,           # skew
    int *,             # i/o info
    int *,
    int *,
    int *,
    int *,
    char *,            # files
    char *,
    char *,
    int *,             # PETSc logging
    int *,
    int *,             # error codes
    char *,
    size_t,            # string lengths
    size_t,
    size_t,
    size_t) except *

cdef extern void id_split "FC_FUNC_(id_split, ID_SPLIT)" (
    int *,
    int *,
    int *,
    int *,
    int *) except *

cdef extern void lnklst "FC_FUNC(lnklst, LNKLST)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void local "FC_FUNC(local, LOCAL)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *) except *

cdef extern void localf "FC_FUNC(localf, LOCALF)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *) except *

cdef extern void localx "FC_FUNC(localx, LOCALX)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *) except *

cdef extern void makemsr "FC_FUNC(makemsr, MAKEMSR)" (
    petsc.Mat *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    double *) except *

cdef extern void matmod_def "FC_FUNC_(matmod_def, MATMOD_DEF)" (int *) except *

cdef extern void nfind "FC_FUNC(nfind, NFIND)" (
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *) except *

cdef extern void preshape "FC_FUNC(preshape, PRESHAPE)" (
    double *,
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void preshape2d "FC_FUNC(preshape2d, PRESHAPE2D)" (
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void read_bc "FC_FUNC_(read_bc, READ_BC)" (
    double *,
    double *,
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void read_connect "FC_FUNC_(read_connect, READ_CONNECT)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *


cdef extern void read_coords "FC_FUNC_(read_coords, READ_COORDS)" (
    double *,
    double *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void read_diff "FC_FUNC_(read_diff, READ_DIFF)" (
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void read_fuldat "FC_FUNC_(read_fuldat, READ_FULDAT)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void read_hist "FC_FUNC_(read_hist, READ_HIST)" (
    double *,
    double *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

# code for reading prestress input files is presently disabled
cdef extern void read_prestr "FC_FUNC_(read_prestr, READ_PRESTR)" ( # MISSING
    double *,
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    char *,
    size_t,
    size_t) except *

cdef extern void read_skew "FC_FUNC_(read_skew, READ_SKEW)" (
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void read_slip "FC_FUNC_(read_slip, READ_SLIP)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void read_split "FC_FUNC_(read_split, READ_SPLIT)" (
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void read_stateout "FC_FUNC_(read_stateout, READ_STATEOUT)" (
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void read_timdat "FC_FUNC_(read_timdat, READ_TIMDAT)" (
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void read_tractions "FC_FUNC_(read_tractions, READ_TRACTIONS)" (
    int *,
    double *,
    double *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void read_wink "FC_FUNC_(read_wink, READ_WINK)" (
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void scan_bc "FC_FUNC_(scan_bc, SCAN_BC)" (
    int *,
    int *,
    char *,
    char *,
    char *,
    char *,
    int *,
    char *,
    size_t,
    size_t,
    size_t,
    size_t,
    size_t) except *

cdef extern void scan_connect "FC_FUNC_(scan_connect, SCAN_CONNECT)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void scan_coords "FC_FUNC_(scan_coords, SCAN_COORDS)" (
    int *,
    int *,
    char *,
    char *,
    int *,
    char *,
    size_t,
    size_t,
    size_t) except *

cdef extern void scan_diff "FC_FUNC_(scan_diff, SCAN_DIFF)" (
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void scan_fuldat "FC_FUNC_(scan_fuldat, SCAN_FULDAT)" (
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void scan_hist "FC_FUNC_(scan_hist, SCAN_HIST)" (
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

# code for reading prestress input files is presently disabled
cdef extern void scan_prestr "FC_FUNC_(scan_prestr, SCAN_PRESTR)" ( # MISSING
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void scan_skew "FC_FUNC_(scan_skew, SCAN_SKEW)" (
    int *,
    int *,
    char *,
    char *,
    int *,
    char *,
    size_t,
    size_t,
    size_t) except *

cdef extern void scan_slip "FC_FUNC_(scan_slip, SCAN_SLIP)" (
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void scan_split "FC_FUNC_(scan_split, SCAN_SPLIT)" (
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void scan_timdat "FC_FUNC_(scan_timdat, SCAN_TIMDAT)" (
    int *,
    int *,
    int *,
    char *,
    char *,
    int *,
    char *,
    size_t,
    size_t,
    size_t) except *

cdef extern void scan_tractions "FC_FUNC_(scan_tractions, SCAN_TRACTIONS)" (
    int *,
    int *,
    int *,
    char *,
    char *,
    int *,
    char *,
    size_t,
    size_t,
    size_t) except *

cdef extern void scan_wink "FC_FUNC_(scan_wink, SCAN_WINK)" (
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void scan_winkx "FC_FUNC_(scan_winkx, SCAN_WINKX)" (
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void sort_elements "FC_FUNC_(sort_elements, SORT_ELEMENTS)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void sort_slip_nodes "FC_FUNC_(sort_slip_nodes, SORT_SLIP_NODES)" (
    int *,
    int *,
    int *,
    int *) except *

cdef extern void sort_split_nodes "FC_FUNC_(sort_split_nodes, SORT_SPLIT_NODES)" (
    int *,
    int *,
    int *,
    int *) except *

cdef extern void try_binio "FC_FUNC_(try_binio, TRY_BINIO)" (
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void viscos "FC_FUNC(viscos, VISCOS)" (
    petsc.Mat *,       # sparse
    petsc.Vec *,
    petsc.Vec *,
    double *,           # force
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,           # 10
    double *,
    double *,
    double *,
    int *,
    double *,
    double *,           # global
    double *,
    double *,
    double *,
    int *,             # 20
    int *,
    double *,
    int *,
    int *,
    int *,             # bc
    double *,
    double *,           # slip
    double *,
    double *,
    double *,           # 30
    int *,
    int *,
    double *,
    int *,
    int *,
    int *,
    double *,           # split
    int *,
    double *,
    double *,           # 40
    double *,           # stiff
    double *,
    double *,           # element
    double *,
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,             # 50
    int *,
    int *,
    int *,
    int *,             # traction
    double *,
    double *,
    double *,
    int *,
    double *,           # material
    int *,             # 60
    double *,           # element type
    double *,
    double *,
    int *,
    double *,           # time data
    double *,
    int *,
    int *,
    int *,
    double *,           # 70
    double *,
    int *,
    int *,
    int *,
    double *,
    double *,
    double *,
    int *,
    double *,           # iterations
    double *,           # 80: skew
    int *,             # i/o info
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,            # files
    char *,
    char *,
    int *,             # 90: PETSC logging
    int *,
    int *,             # error codes
    char *,
    size_t,            # string lengths
    size_t,
    size_t,
    size_t) except *

cdef extern void viscos_setup "FC_FUNC_(viscos_setup, VISCOS_SETUP)" (
    int *,
    int *,
    char *,
    char *,
    int *,
    int *,
    char *,
    size_t,
    size_t,
    size_t) except *

cdef extern void viscos_step "FC_FUNC_(viscos_step, VISCOS_STEP)" (
    petsc.Mat *,       # sparse
    petsc.Vec *,
    petsc.Vec *,
    double *,           # force
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,           # 10
    double *,
    double *,
    double *,
    int *,
    double *,
    double *,           # global
    double *,
    double *,
    double *,
    int *,             # 20
    int *,
    double *,
    int *,
    int *,
    int *,             # bc
    double *,
    double *,           # slip
    double *,
    double *,
    double *,           # 30
    int *,
    int *,
    double *,
    int *,
    int *,
    int *,
    double *,           # split
    int *,
    double *,
    double *,           # 40
    double *,           # stiff
    double *,
    double *,           # element
    double *,
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,             # 50
    int *,
    int *,
    int *,
    int *,             # traction
    double *,
    double *,
    double *,
    int *,
    double *,           # material
    int *,             # 60
    double *,           # element type
    double *,
    double *,
    int *,
    double *,           # time data
    double *,
    int *,
    int *,
    int *,
    double *,           # 70
    double *,
    int *,
    int *,
    int *,
    double *,
    double *,
    double *,
    int *,
    double *,           # iterations
    double *,           # 80: skew
    int *,             # i/o info
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,            # files
    char *,
    char *,
    int *,             # 90: PETSC logging
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    double *,
    double *,
    int *,
    double *,
    int *,             # error codes
    char *,
    size_t,            # string lengths
    size_t,
    size_t,
    size_t) except *

cdef extern void viscos_cleanup "FC_FUNC_(viscos_cleanup, VISCOS_CLEANUP)" (
    int *,
    int *,
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void write_bc "FC_FUNC_(write_bc, WRITE_BC)" (
    double *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void write_connect "FC_FUNC_(write_connect, WRITE_CONNECT)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    char *,
    int *,
    char *,
    size_t,
    size_t,
    size_t) except *

cdef extern void write_coords "FC_FUNC_(write_coords, WRITE_COORDS)" (
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    char *,
    int *,
    char *,
    size_t,
    size_t,
    size_t) except *

cdef extern void write_diff "FC_FUNC_(write_diff, WRITE_DIFF)" (
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void write_element_info "FC_FUNC_(write_element_info, WRITE_ELEMENT_INFO)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    double *,
    double *,
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void write_fuldat "FC_FUNC_(write_fuldat, WRITE_FULDAT)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    char *,
    int *,
    char *,
    size_t,
    size_t,
    size_t) except *

cdef extern void write_global_info "FC_FUNC_(write_global_info, WRITE_GLOBAL_INFO)" (
    char *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    char *,
    size_t,
    size_t,
    size_t) except *

cdef extern void write_hist "FC_FUNC_(write_hist, WRITE_HIST)" (
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void write_props "FC_FUNC_(write_props, WRITE_PROPS)" (
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    char *,
    int *,
    char *,
    size_t,
    size_t,
    size_t) except *

cdef extern void write_skew "FC_FUNC_(write_skew, WRITE_SKEW)" (
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void write_slip "FC_FUNC_(write_slip, WRITE_SLIP)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    char *,
    int *,
    char *,
    size_t,
    size_t,
    size_t) except *

cdef extern void write_sparse_info "FC_FUNC_(write_sparse_info, WRITE_SPARSE_INFO)" (
    int *,
    int *,
    int *,
    int *,
    double *,
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void write_split "FC_FUNC_(write_split, WRITE_SPLIT)" (
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    char *,
    int *,
    char *,
    size_t,
    size_t,
    size_t) except *

cdef extern void write_split_plot "FC_FUNC_(write_split_plot, WRITE_SPLIT_PLOT)" (
    int *,
    int *,
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void write_stateout "FC_FUNC_(write_stateout, WRITE_STATEOUT)" (
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    char *,
    int *,
    char *,
    size_t,
    size_t,
    size_t) except *

cdef extern void write_strscomp "FC_FUNC_(write_strscomp, WRITE_STRSCOMP)" (
    double *,
    double *,
    double *,
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void write_subiter "FC_FUNC_(write_subiter, WRITE_SUBITER)" (
    int *,
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void write_tractions "FC_FUNC_(write_tractions, WRITE_TRACTIONS)" (
    int *,
    double *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void write_timdat "FC_FUNC_(write_timdat, WRITE_TIMDAT)" (
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void write_ucd_mesh "FC_FUNC_(write_ucd_mesh, WRITE_UCD_MESH)" (
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    size_t) except *

cdef extern void write_wink "FC_FUNC_(write_wink, WRITE_WINK)" (
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *

cdef extern void write_winkx "FC_FUNC_(write_winkx, WRITE_WINKX)" (
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    char *,
    int *,
    char *,
    size_t,
    size_t) except *


# end of file
