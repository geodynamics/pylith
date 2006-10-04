// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
//
//  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
//
//  Permission is hereby granted, free of charge, to any person obtaining
//  a copy of this software and associated documentation files (the
//  "Software"), to deal in the Software without restriction, including
//  without limitation the rights to use, copy, modify, merge, publish,
//  distribute, sublicense, and/or sell copies of the Software, and to
//  permit persons to whom the Software is furnished to do so, subject to
//  the following conditions:
//
//  The above copyright notice and this permission notice shall be
//  included in all copies or substantial portions of the Software.
//
//  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
//  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
//  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
//  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
//  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
//  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
#include "sorting.h"        // element sorting routines
#include "viscos.h"        // time-dependent solution driver
#include "write_modelinfo.h" // output routines
#include "misc.h"          // miscellaneous methods

// the method table

struct PyMethodDef pypylith3d_methods[] = {

    // initialize PETSc
    {pypylith3d_petsc_initialize__name__, pypylith3d_petsc_initialize,
     METH_VARARGS, pypylith3d_petsc_initialize__doc__},

    // finalize PETSc
    {pypylith3d_petsc_finalize__name__, pypylith3d_petsc_finalize,
     METH_VARARGS, pypylith3d_petsc_finalize__doc__},

    // Setup PETSc Logging
    {pypylith3d_setup_petsc_logging__name__, pypylith3d_setup_petsc_logging,
     METH_VARARGS, pypylith3d_setup_petsc_logging__doc__},

    // allocate an integer array
    {pypylith3d_allocateInt__name__, pypylith3d_allocateInt,
     METH_VARARGS, pypylith3d_allocateInt__doc__},

    // allocate a double array
    {pypylith3d_allocateDouble__name__, pypylith3d_allocateDouble,
     METH_VARARGS, pypylith3d_allocateDouble__doc__},

    // assign equation numbers for Winkler forces
    {pypylith3d_assign_wink__name__, pypylith3d_assign_wink,
     METH_VARARGS, pypylith3d_assign_wink__doc__},

    // compute gravitational prestresses
    {pypylith3d_autoprestr__name__, pypylith3d_autoprestr,
     METH_VARARGS, pypylith3d_autoprestr__doc__},

    // compute maximum number of nonzero entries in stiffness matrix
    {pypylith3d_cmp_stiffsz__name__, pypylith3d_cmp_stiffsz,
     METH_VARARGS, pypylith3d_cmp_stiffsz__doc__},

    // create id and idx arrays
    {pypylith3d_create_id__name__, pypylith3d_create_id,
     METH_VARARGS, pypylith3d_create_id__doc__},

    // convert a double list to an array
    {pypylith3d_doubleListToArray__name__, pypylith3d_doubleListToArray,
     METH_VARARGS, pypylith3d_doubleListToArray__doc__},

    // Retrieve an integer list member
    {pypylith3d_intListRef__name__, pypylith3d_intListRef,
     METH_VARARGS, pypylith3d_intListRef__doc__},

    // Retrieve a double list member
    {pypylith3d_doubleListRef__name__, pypylith3d_doubleListRef,
     METH_VARARGS, pypylith3d_doubleListRef__doc__},

    // Set an integer list member
    {pypylith3d_intListSet__name__, pypylith3d_intListSet,
     METH_VARARGS, pypylith3d_intListSet__doc__},

    // Set a double list member
    {pypylith3d_doubleListSet__name__, pypylith3d_doubleListSet,
     METH_VARARGS, pypylith3d_doubleListSet__doc__},

    // drive elastic solution
    {pypylith3d_elastc__name__, pypylith3d_elastc,
     METH_VARARGS, pypylith3d_elastc__doc__},

    // form id array for split nodes
    {pypylith3d_id_split__name__, pypylith3d_id_split,
     METH_VARARGS, pypylith3d_id_split__doc__},

    // convert an integer list to an array
    {pypylith3d_intListToArray__name__, pypylith3d_intListToArray,
     METH_VARARGS, pypylith3d_intListToArray__doc__},

    // create linked list for sparse matrix
    {pypylith3d_lnklst__name__, pypylith3d_lnklst,
     METH_VARARGS, pypylith3d_lnklst__doc__},

    // localize id array
    {pypylith3d_local__name__, pypylith3d_local,
     METH_VARARGS, pypylith3d_local__doc__},

    // localize id array for split nodes
    {pypylith3d_localf__name__, pypylith3d_localf,
     METH_VARARGS, pypylith3d_localf__doc__},

    // localize id array for slippery nodes
    {pypylith3d_localx__name__, pypylith3d_localx,
     METH_VARARGS, pypylith3d_localx__doc__},

    // create PETSc matrix
    {pypylith3d_createPETScMat__name__, pypylith3d_createPETScMat,
     METH_VARARGS, pypylith3d_createPETScMat__doc__},

    // destroy PETSc matrix
    {pypylith3d_destroyPETScMat__name__, pypylith3d_destroyPETScMat,
     METH_VARARGS, pypylith3d_destroyPETScMat__doc__},

    // output PETSc Mesh and Fields
    {pypylith3d_outputMesh__name__, pypylith3d_outputMesh,
     METH_VARARGS, pypylith3d_outputMesh__doc__},

    // create sparse matrix in modified sparse row format
    {pypylith3d_makemsr__name__, pypylith3d_makemsr,
     METH_VARARGS, pypylith3d_makemsr__doc__},

    // initialize material model information
    {pypylith3d_matmod_def__name__, pypylith3d_matmod_def,
     METH_VARARGS, pypylith3d_matmod_def__doc__},

    // find closest fault neignbors for slippery nodes
    {pypylith3d_nfind__name__, pypylith3d_nfind,
     METH_VARARGS, pypylith3d_nfind__doc__},

    // precompute shape function information
    {pypylith3d_preshape__name__, pypylith3d_preshape,
     METH_VARARGS, pypylith3d_preshape__doc__},

    // precompute shape function information for surfaces
    {pypylith3d_preshape2d__name__, pypylith3d_preshape2d,
     METH_VARARGS, pypylith3d_preshape2d__doc__},

    // read boundary conditions
    {pypylith3d_read_bc__name__, pypylith3d_read_bc,
     METH_VARARGS, pypylith3d_read_bc__doc__},

    // read connectivities
    {pypylith3d_read_connect__name__, pypylith3d_read_connect,
     METH_VARARGS, pypylith3d_read_connect__doc__},

    // read coordinates
    {pypylith3d_read_coords__name__, pypylith3d_read_coords,
     METH_VARARGS, pypylith3d_read_coords__doc__},

    // read differential forces
    {pypylith3d_read_diff__name__, pypylith3d_read_diff,
     METH_VARARGS, pypylith3d_read_diff__doc__},

    // read time steps at which full output is desired
    {pypylith3d_read_fuldat__name__, pypylith3d_read_fuldat,
     METH_VARARGS, pypylith3d_read_fuldat__doc__},

    // read time history info
    {pypylith3d_read_hist__name__, pypylith3d_read_hist,
     METH_VARARGS, pypylith3d_read_hist__doc__},

    // read prestresses
    // {pypylith3d_read_prestr__name__, pypylith3d_read_prestr,
     // METH_VARARGS, pypylith3d_read_prestr__doc__},

    // read local coordinate rotations
    {pypylith3d_read_skew__name__, pypylith3d_read_skew,
     METH_VARARGS, pypylith3d_read_skew__doc__},

    // read slippery node definitions
    {pypylith3d_read_slip__name__, pypylith3d_read_slip,
     METH_VARARGS, pypylith3d_read_slip__doc__},

    // read split node definitions
    {pypylith3d_read_split__name__, pypylith3d_read_split,
     METH_VARARGS, pypylith3d_read_split__doc__},

    // read state output information
    {pypylith3d_read_stateout__name__, pypylith3d_read_stateout,
     METH_VARARGS, pypylith3d_read_stateout__doc__},

    // read time step group information
    {pypylith3d_read_timdat__name__, pypylith3d_read_timdat,
     METH_VARARGS, pypylith3d_read_timdat__doc__},

    // read traction boundary conditions
    // {pypylith3d_read_traction__name__, pypylith3d_read_traction,
     // METH_VARARGS, pypylith3d_read_traction__doc__},

    // read winkler force information
    {pypylith3d_read_wink__name__, pypylith3d_read_wink,
     METH_VARARGS, pypylith3d_read_wink__doc__},

    // process mesh
    {pypylith3d_processMesh__name__, pypylith3d_processMesh,
     METH_VARARGS, pypylith3d_processMesh__doc__},

    // scan boundary condition file
    {pypylith3d_scan_bc__name__, pypylith3d_scan_bc,
     METH_VARARGS, pypylith3d_scan_bc__doc__},

    // scan connectivity file
    {pypylith3d_scan_connect__name__, pypylith3d_scan_connect,
     METH_VARARGS, pypylith3d_scan_connect__doc__},

    // scan coordinates file
    {pypylith3d_scan_coords__name__, pypylith3d_scan_coords,
     METH_VARARGS, pypylith3d_scan_coords__doc__},

    // scan differential forces file
    {pypylith3d_scan_diff__name__, pypylith3d_scan_diff,
     METH_VARARGS, pypylith3d_scan_diff__doc__},

    // scan time step output info file
    {pypylith3d_scan_fuldat__name__, pypylith3d_scan_fuldat,
     METH_VARARGS, pypylith3d_scan_fuldat__doc__},

    // scan time history definition file
    {pypylith3d_scan_hist__name__, pypylith3d_scan_hist,
     METH_VARARGS, pypylith3d_scan_hist__doc__},

    // scan prestress file
    // {pypylith3d_scan_prestr__name__, pypylith3d_scan_prestr,
     // METH_VARARGS, pypylith3d_scan_prestr__doc__},

    // scan local coordinate rotations file
    {pypylith3d_scan_skew__name__, pypylith3d_scan_skew,
     METH_VARARGS, pypylith3d_scan_skew__doc__},

    // scan slippery node definitions file
    {pypylith3d_scan_slip__name__, pypylith3d_scan_slip,
     METH_VARARGS, pypylith3d_scan_slip__doc__},

    // scan split node definitions file
    {pypylith3d_scan_split__name__, pypylith3d_scan_split,
     METH_VARARGS, pypylith3d_scan_split__doc__},

    // scan time step group info file
    {pypylith3d_scan_timdat__name__, pypylith3d_scan_timdat,
     METH_VARARGS, pypylith3d_scan_timdat__doc__},

    // scan traction boundary conditions file
    // {pypylith3d_scan_tractions__name__, pypylith3d_scan_tractions,
     // METH_VARARGS, pypylith3d_scan_tractions__doc__},

    // scan winkler forces info file
    {pypylith3d_scan_wink__name__, pypylith3d_scan_wink,
     METH_VARARGS, pypylith3d_scan_wink__doc__},

    // scan winkler forces info file for slippery nodes
    {pypylith3d_scan_winkx__name__, pypylith3d_scan_winkx,
     METH_VARARGS, pypylith3d_scan_winkx__doc__},

    // sort elements into families
    {pypylith3d_sort_elements__name__, pypylith3d_sort_elements,
     METH_VARARGS, pypylith3d_sort_elements__doc__},

    // sort slippery node elements
    {pypylith3d_sort_slip_nodes__name__, pypylith3d_sort_slip_nodes,
     METH_VARARGS, pypylith3d_sort_slip_nodes__doc__},

    // sort split node elements
    {pypylith3d_sort_split_nodes__name__, pypylith3d_sort_split_nodes,
     METH_VARARGS, pypylith3d_sort_split_nodes__doc__},

    // drive time-dependent solution
    {pypylith3d_viscos__name__, pypylith3d_viscos,
     METH_VARARGS, pypylith3d_viscos__doc__},

    // setup time-dependent solution
    {pypylith3d_viscos_setup__name__, pypylith3d_viscos_setup,
     METH_VARARGS, pypylith3d_viscos_setup__doc__},

    // drive time-dependent solution
    {pypylith3d_viscos_step__name__, pypylith3d_viscos_step,
     METH_VARARGS, pypylith3d_viscos_step__doc__},

    // cleanup time-dependent solution
    {pypylith3d_viscos_cleanup__name__, pypylith3d_viscos_cleanup,
     METH_VARARGS, pypylith3d_viscos_cleanup__doc__},

    // write out BC
    {pypylith3d_write_bc__name__, pypylith3d_write_bc,
     METH_VARARGS, pypylith3d_write_bc__doc__},

    // write out element node array
    {pypylith3d_write_connect__name__, pypylith3d_write_connect,
     METH_VARARGS, pypylith3d_write_connect__doc__},

    // write out nodal coordinates
    {pypylith3d_write_coords__name__, pypylith3d_write_coords,
     METH_VARARGS, pypylith3d_write_coords__doc__},

    // write out differential forces
    {pypylith3d_write_diff__name__, pypylith3d_write_diff,
     METH_VARARGS, pypylith3d_write_diff__doc__},

    // write out element and prestress information
    {pypylith3d_write_element_info__name__, pypylith3d_write_element_info,
     METH_VARARGS, pypylith3d_write_element_info__doc__},

    // write out time steps at which full output is desired
    {pypylith3d_write_fuldat__name__, pypylith3d_write_fuldat,
     METH_VARARGS, pypylith3d_write_fuldat__doc__},

    // write out global control information
    {pypylith3d_write_global_info__name__, pypylith3d_write_global_info,
     METH_VARARGS, pypylith3d_write_global_info__doc__},

    // write out time histories
    {pypylith3d_write_hist__name__, pypylith3d_write_hist,
     METH_VARARGS, pypylith3d_write_hist__doc__},

    // write out material property information
    {pypylith3d_write_props__name__, pypylith3d_write_props,
     METH_VARARGS, pypylith3d_write_props__doc__},

    // write out local coordinate rotations
    {pypylith3d_write_skew__name__, pypylith3d_write_skew,
     METH_VARARGS, pypylith3d_write_skew__doc__},

    // write out slippery node definitions
    {pypylith3d_write_slip__name__, pypylith3d_write_slip,
     METH_VARARGS, pypylith3d_write_slip__doc__},

    // write out sparse matrix information
    {pypylith3d_write_sparse_info__name__, pypylith3d_write_sparse_info,
     METH_VARARGS, pypylith3d_write_sparse_info__doc__},

    // write out split node definitions
    {pypylith3d_write_split__name__, pypylith3d_write_split,
     METH_VARARGS, pypylith3d_write_split__doc__},

    // write out split node info for plotting
    {pypylith3d_write_split_plot__name__, pypylith3d_write_split_plot,
     METH_VARARGS, pypylith3d_write_split_plot__doc__},

    // write out state variable output info
    {pypylith3d_write_stateout__name__, pypylith3d_write_stateout,
     METH_VARARGS, pypylith3d_write_stateout__doc__},

    // write out stress computation information
    {pypylith3d_write_strscomp__name__, pypylith3d_write_strscomp,
     METH_VARARGS, pypylith3d_write_strscomp__doc__},

    // write out subiteration convergence information
    {pypylith3d_write_subiter__name__, pypylith3d_write_subiter,
     METH_VARARGS, pypylith3d_write_subiter__doc__},

    // write out time step information
    {pypylith3d_write_timdat__name__, pypylith3d_write_timdat,
     METH_VARARGS, pypylith3d_write_timdat__doc__},

    // write out traction information
    {pypylith3d_write_tractions__name__, pypylith3d_write_tractions,
     METH_VARARGS, pypylith3d_write_tractions__doc__},

    // write mesh info to UCD file
    {pypylith3d_write_ucd_mesh__name__, pypylith3d_write_ucd_mesh,
     METH_VARARGS, pypylith3d_write_ucd_mesh__doc__},

    // write out Winkler force info
    {pypylith3d_write_wink__name__, pypylith3d_write_wink,
     METH_VARARGS, pypylith3d_write_wink__doc__},

    // write out slippery node Winkler force info
    {pypylith3d_write_winkx__name__, pypylith3d_write_winkx,
     METH_VARARGS, pypylith3d_write_winkx__doc__},

    // try binary I/O to see if it works
    {pypylith3d_try_binio__name__, pypylith3d_try_binio,
     METH_VARARGS, pypylith3d_try_binio__doc__},

    // copyright note
    {pypylith3d_copyright__name__, pypylith3d_copyright,
     METH_VARARGS, pypylith3d_copyright__doc__},


// Sentinel
    {0, 0}
};

// version
// $Id: bindings.cc,v 1.10 2005/04/21 23:16:35 willic3 Exp $

// End of file
