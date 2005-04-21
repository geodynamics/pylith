// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                               Charles A. Williams
//                        Rensselaer Polytechnic Institute
//                        (C) 2005 All Rights Reserved
// 
//  All worldwide rights reserved.  A license to use, copy, modify and
//  distribute this software for non-commercial research purposes only
//  is hereby granted, provided that this copyright notice and
//  accompanying disclaimer is not modified or removed from the software.
//
//  DISCLAIMER:  The software is distributed "AS IS" without any express
//  or implied warranty, including but not limited to, any implied
//  warranties of merchantability or fitness for a particular purpose
//  or any warranty of non-infringement of any current or pending patent
//  rights.  The authors of the software make no representations about
//  the suitability of this software for any particular purpose.  The
//  entire risk as to the quality and performance of the software is with
//  the user.  Should the software prove defective, the user assumes the
//  cost of all necessary servicing, repair or correction.  In
//  particular, neither Rensselaer Polytechnic Institute, nor the authors
//  of the software are liable for any indirect, special, consequential,
//  or incidental damages related to the software, to the maximum extent
//  the law permits.
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 

#include <portinfo>
#include <Python.h>

#include "bindings.h"

#include "array.h"         // array allocation and conversion functions
#include "autoprestr.h"    // prestress autocomputation
#include "elastc.h"        // elastic solution driver
#include "numbering.h"     // routines to create global equation numbers
                           // and localize them.
#include "parser.h"        // parsers
#include "scanner.h"       // scanners
#include "setup.h"         // initialization/setup routines
#include "sparse.h"        // sparse matrix routines
#include "viscos.h"        // time-dependent solution driver
#include "write_modelinfo.h" // output routines
#include "misc.h"          // miscellaneous methods

// the method table

struct PyMethodDef pylithomop3d_methods[] = {

    // initialize PETSc
    {pylithomop3d_petsc_initialize__name__, pylithomop3d_petsc_initialize,
     METH_VARARGS, pylithomop3d_petsc_initialize__doc__},

    // finalize PETSc
    {pylithomop3d_petsc_finalize__name__, pylithomop3d_petsc_finalize,
     METH_VARARGS, pylithomop3d_petsc_finalize__doc__},

    // Setup PETSc Logging
    {pylithomop3d_setup_petsc_logging__name__, pylithomop3d_setup_petsc_logging,
     METH_VARARGS, pylithomop3d_setup_petsc_logging__doc__},

    // allocate an integer array
    {pylithomop3d_allocateInt__name__, pylithomop3d_allocateInt,
     METH_VARARGS, pylithomop3d_allocateInt__doc__},

    // allocate a double array
    {pylithomop3d_allocateDouble__name__, pylithomop3d_allocateDouble,
     METH_VARARGS, pylithomop3d_allocateDouble__doc__},

    // assign equation numbers for Winkler forces
    {pylithomop3d_assign_wink__name__, pylithomop3d_assign_wink,
     METH_VARARGS, pylithomop3d_assign_wink__doc__},

    // compute gravitational prestresses
    {pylithomop3d_autoprestr__name__, pylithomop3d_autoprestr,
     METH_VARARGS, pylithomop3d_autoprestr__doc__},

    // compute maximum number of nonzero entries in stiffness matrix
    {pylithomop3d_cmp_stiffsz__name__, pylithomop3d_cmp_stiffsz,
     METH_VARARGS, pylithomop3d_cmp_stiffsz__doc__},

    // create id and idx arrays
    {pylithomop3d_create_id__name__, pylithomop3d_create_id,
     METH_VARARGS, pylithomop3d_create_id__doc__},

    // convert a double list to an array
    {pylithomop3d_doubleListToArray__name__, pylithomop3d_doubleListToArray,
     METH_VARARGS, pylithomop3d_doubleListToArray__doc__},

    // drive elastic solution
    {pylithomop3d_elastc__name__, pylithomop3d_elastc,
     METH_VARARGS, pylithomop3d_elastc__doc__},

    // form id array for split nodes
    {pylithomop3d_id_split__name__, pylithomop3d_id_split,
     METH_VARARGS, pylithomop3d_id_split__doc__},

    // convert an integer list to an array
    {pylithomop3d_intListToArray__name__, pylithomop3d_intListToArray,
     METH_VARARGS, pylithomop3d_intListToArray__doc__},

    // create linked list for sparse matrix
    {pylithomop3d_lnklst__name__, pylithomop3d_lnklst,
     METH_VARARGS, pylithomop3d_lnklst__doc__},

    // localize id array
    {pylithomop3d_local__name__, pylithomop3d_local,
     METH_VARARGS, pylithomop3d_local__doc__},

    // localize id array for split nodes
    {pylithomop3d_localf__name__, pylithomop3d_localf,
     METH_VARARGS, pylithomop3d_localf__doc__},

    // localize id array for slippery nodes
    {pylithomop3d_localx__name__, pylithomop3d_localx,
     METH_VARARGS, pylithomop3d_localx__doc__},

    // create PETSc matrix
    {pylithomop3d_createPETScMat__name__, pylithomop3d_createPETScMat,
     METH_VARARGS, pylithomop3d_createPETScMat__doc__},

    // destroy PETSc matrix
    {pylithomop3d_destroyPETScMat__name__, pylithomop3d_destroyPETScMat,
     METH_VARARGS, pylithomop3d_destroyPETScMat__doc__},

    // create sparse matrix in modified sparse row format
    {pylithomop3d_makemsr__name__, pylithomop3d_makemsr,
     METH_VARARGS, pylithomop3d_makemsr__doc__},

    // initialize material model information
    {pylithomop3d_matmod_def__name__, pylithomop3d_matmod_def,
     METH_VARARGS, pylithomop3d_matmod_def__doc__},

    // find closest fault neignbors for slippery nodes
    {pylithomop3d_nfind__name__, pylithomop3d_nfind,
     METH_VARARGS, pylithomop3d_nfind__doc__},

    // precompute shape function information
    {pylithomop3d_preshape__name__, pylithomop3d_preshape,
     METH_VARARGS, pylithomop3d_preshape__doc__},

    // read boundary conditions
    {pylithomop3d_read_bc__name__, pylithomop3d_read_bc,
     METH_VARARGS, pylithomop3d_read_bc__doc__},

    // read connectivities
    {pylithomop3d_read_connect__name__, pylithomop3d_read_connect,
     METH_VARARGS, pylithomop3d_read_connect__doc__},

    // read coordinates
    {pylithomop3d_read_coords__name__, pylithomop3d_read_coords,
     METH_VARARGS, pylithomop3d_read_coords__doc__},

    // read differential forces
    {pylithomop3d_read_diff__name__, pylithomop3d_read_diff,
     METH_VARARGS, pylithomop3d_read_diff__doc__},

    // read time steps at which full output is desired
    {pylithomop3d_read_fuldat__name__, pylithomop3d_read_fuldat,
     METH_VARARGS, pylithomop3d_read_fuldat__doc__},

    // read time history info
    {pylithomop3d_read_hist__name__, pylithomop3d_read_hist,
     METH_VARARGS, pylithomop3d_read_hist__doc__},

    // read prestresses
    // {pylithomop3d_read_prestr__name__, pylithomop3d_read_prestr,
     // METH_VARARGS, pylithomop3d_read_prestr__doc__},

    // read local coordinate rotations
    {pylithomop3d_read_skew__name__, pylithomop3d_read_skew,
     METH_VARARGS, pylithomop3d_read_skew__doc__},

    // read slippery node definitions
    {pylithomop3d_read_slip__name__, pylithomop3d_read_slip,
     METH_VARARGS, pylithomop3d_read_slip__doc__},

    // read split node definitions
    {pylithomop3d_read_split__name__, pylithomop3d_read_split,
     METH_VARARGS, pylithomop3d_read_split__doc__},

    // read state output information
    {pylithomop3d_read_stateout__name__, pylithomop3d_read_stateout,
     METH_VARARGS, pylithomop3d_read_stateout__doc__},

    // read time step group information
    {pylithomop3d_read_timdat__name__, pylithomop3d_read_timdat,
     METH_VARARGS, pylithomop3d_read_timdat__doc__},

    // read traction boundary conditions
    // {pylithomop3d_read_traction__name__, pylithomop3d_read_traction,
     // METH_VARARGS, pylithomop3d_read_traction__doc__},

    // read winkler force information
    {pylithomop3d_read_wink__name__, pylithomop3d_read_wink,
     METH_VARARGS, pylithomop3d_read_wink__doc__},

    // scan boundary condition file
    {pylithomop3d_scan_bc__name__, pylithomop3d_scan_bc,
     METH_VARARGS, pylithomop3d_scan_bc__doc__},

    // scan connectivity file
    {pylithomop3d_scan_connect__name__, pylithomop3d_scan_connect,
     METH_VARARGS, pylithomop3d_scan_connect__doc__},

    // scan coordinates file
    {pylithomop3d_scan_coords__name__, pylithomop3d_scan_coords,
     METH_VARARGS, pylithomop3d_scan_coords__doc__},

    // scan differential forces file
    {pylithomop3d_scan_diff__name__, pylithomop3d_scan_diff,
     METH_VARARGS, pylithomop3d_scan_diff__doc__},

    // scan time step output info file
    {pylithomop3d_scan_fuldat__name__, pylithomop3d_scan_fuldat,
     METH_VARARGS, pylithomop3d_scan_fuldat__doc__},

    // scan time history definition file
    {pylithomop3d_scan_hist__name__, pylithomop3d_scan_hist,
     METH_VARARGS, pylithomop3d_scan_hist__doc__},

    // scan prestress file
    // {pylithomop3d_scan_prestr__name__, pylithomop3d_scan_prestr,
     // METH_VARARGS, pylithomop3d_scan_prestr__doc__},

    // scan local coordinate rotations file
    {pylithomop3d_scan_skew__name__, pylithomop3d_scan_skew,
     METH_VARARGS, pylithomop3d_scan_skew__doc__},

    // scan slippery node definitions file
    {pylithomop3d_scan_slip__name__, pylithomop3d_scan_slip,
     METH_VARARGS, pylithomop3d_scan_slip__doc__},

    // scan split node definitions file
    {pylithomop3d_scan_split__name__, pylithomop3d_scan_split,
     METH_VARARGS, pylithomop3d_scan_split__doc__},

    // scan time step group info file
    {pylithomop3d_scan_timdat__name__, pylithomop3d_scan_timdat,
     METH_VARARGS, pylithomop3d_scan_timdat__doc__},

    // scan traction boundary conditions file
    // {pylithomop3d_scan_traction__name__, pylithomop3d_scan_traction,
     // METH_VARARGS, pylithomop3d_scan_traction__doc__},

    // scan winkler forces info file
    {pylithomop3d_scan_wink__name__, pylithomop3d_scan_wink,
     METH_VARARGS, pylithomop3d_scan_wink__doc__},

    // scan winkler forces info file for slippery nodes
    {pylithomop3d_scan_winkx__name__, pylithomop3d_scan_winkx,
     METH_VARARGS, pylithomop3d_scan_winkx__doc__},

    // sort elements into families
    {pylithomop3d_sort_elements__name__, pylithomop3d_sort_elements,
     METH_VARARGS, pylithomop3d_sort_elements__doc__},

    // sort slippery node elements
    {pylithomop3d_sort_slip_nodes__name__, pylithomop3d_sort_slip_nodes,
     METH_VARARGS, pylithomop3d_sort_slip_nodes__doc__},

    // sort split node elements
    {pylithomop3d_sort_split_nodes__name__, pylithomop3d_sort_split_nodes,
     METH_VARARGS, pylithomop3d_sort_split_nodes__doc__},

    // drive time-dependent solution
    {pylithomop3d_viscos__name__, pylithomop3d_viscos,
     METH_VARARGS, pylithomop3d_viscos__doc__},

    // write out BC
    {pylithomop3d_write_bc__name__, pylithomop3d_write_bc,
     METH_VARARGS, pylithomop3d_write_bc__doc__},

    // write out element node array
    {pylithomop3d_write_connect__name__, pylithomop3d_write_connect,
     METH_VARARGS, pylithomop3d_write_connect__doc__},

    // write out nodal coordinates
    {pylithomop3d_write_coords__name__, pylithomop3d_write_coords,
     METH_VARARGS, pylithomop3d_write_coords__doc__},

    // write out differential forces
    {pylithomop3d_write_diff__name__, pylithomop3d_write_diff,
     METH_VARARGS, pylithomop3d_write_diff__doc__},

    // write out element and prestress information
    {pylithomop3d_write_element_info__name__, pylithomop3d_write_element_info,
     METH_VARARGS, pylithomop3d_write_element_info__doc__},

    // write out time steps at which full output is desired
    {pylithomop3d_write_fuldat__name__, pylithomop3d_write_fuldat,
     METH_VARARGS, pylithomop3d_write_fuldat__doc__},

    // write out global control information
    {pylithomop3d_write_global_info__name__, pylithomop3d_write_global_info,
     METH_VARARGS, pylithomop3d_write_global_info__doc__},

    // write out time histories
    {pylithomop3d_write_hist__name__, pylithomop3d_write_hist,
     METH_VARARGS, pylithomop3d_write_hist__doc__},

    // write out material property information
    {pylithomop3d_write_props__name__, pylithomop3d_write_props,
     METH_VARARGS, pylithomop3d_write_props__doc__},

    // write out local coordinate rotations
    {pylithomop3d_write_skew__name__, pylithomop3d_write_skew,
     METH_VARARGS, pylithomop3d_write_skew__doc__},

    // write out slippery node definitions
    {pylithomop3d_write_slip__name__, pylithomop3d_write_slip,
     METH_VARARGS, pylithomop3d_write_slip__doc__},

    // write out sparse matrix information
    {pylithomop3d_write_sparse_info__name__, pylithomop3d_write_sparse_info,
     METH_VARARGS, pylithomop3d_write_sparse_info__doc__},

    // write out split node definitions
    {pylithomop3d_write_split__name__, pylithomop3d_write_split,
     METH_VARARGS, pylithomop3d_write_split__doc__},

    // write out split node info for plotting
    {pylithomop3d_write_split_plot__name__, pylithomop3d_write_split_plot,
     METH_VARARGS, pylithomop3d_write_split_plot__doc__},

    // write out state variable output info
    {pylithomop3d_write_stateout__name__, pylithomop3d_write_stateout,
     METH_VARARGS, pylithomop3d_write_stateout__doc__},

    // write out stress computation information
    {pylithomop3d_write_strscomp__name__, pylithomop3d_write_strscomp,
     METH_VARARGS, pylithomop3d_write_strscomp__doc__},

    // write out subiteration convergence information
    {pylithomop3d_write_subiter__name__, pylithomop3d_write_subiter,
     METH_VARARGS, pylithomop3d_write_subiter__doc__},

    // write out time step information
    {pylithomop3d_write_timdat__name__, pylithomop3d_write_timdat,
     METH_VARARGS, pylithomop3d_write_timdat__doc__},

    // write mesh info to UCD file
    {pylithomop3d_write_ucd_mesh__name__, pylithomop3d_write_ucd_mesh,
     METH_VARARGS, pylithomop3d_write_ucd_mesh__doc__},

    // write out Winkler force info
    {pylithomop3d_write_wink__name__, pylithomop3d_write_wink,
     METH_VARARGS, pylithomop3d_write_wink__doc__},

    // write out slippery node Winkler force info
    {pylithomop3d_write_winkx__name__, pylithomop3d_write_winkx,
     METH_VARARGS, pylithomop3d_write_winkx__doc__},


    // copyright note
    {pylithomop3d_copyright__name__, pylithomop3d_copyright,
     METH_VARARGS, pylithomop3d_copyright__doc__},


// Sentinel
    {0, 0}
};

// version
// $Id: bindings.cc,v 1.9 2005/04/21 01:00:30 willic3 Exp $

// End of file
