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

#if !defined(pypylith3d_libpylith3d_h)
#define pypylith3d_libpylith3d_h


// Compute gravitational prestresses.
extern char pypylith3d_autoprestr__name__[];
extern char pypylith3d_autoprestr__doc__[];
extern "C"
PyObject * pypylith3d_autoprestr(PyObject *, PyObject *);


// drive elastic solution
extern char pypylith3d_elastc__name__[];
extern char pypylith3d_elastc__doc__[];
extern "C"
PyObject * pypylith3d_elastc(PyObject *, PyObject *);


// assign equation numbers for Winkler BC
extern char pypylith3d_assign_wink__name__[];
extern char pypylith3d_assign_wink__doc__[];
extern "C"
PyObject * pypylith3d_assign_wink(PyObject *, PyObject *);

// create id array
extern char pypylith3d_create_id__name__[];
extern char pypylith3d_create_id__doc__[];
extern "C"
PyObject * pypylith3d_create_id(PyObject *, PyObject *);

// form id array for split nodes
extern char pypylith3d_id_split__name__[];
extern char pypylith3d_id_split__doc__[];
extern "C"
PyObject * pypylith3d_id_split(PyObject *, PyObject *);

// find closest fault neignbors for slippery nodes
extern char pypylith3d_nfind__name__[];
extern char pypylith3d_nfind__doc__[];
extern "C"
PyObject * pypylith3d_nfind(PyObject *, PyObject *);


// read boundary conditions
extern char pypylith3d_read_bc__name__[];
extern char pypylith3d_read_bc__doc__[];
extern "C"
PyObject * pypylith3d_read_bc(PyObject *, PyObject *);

// read connectivities
extern char pypylith3d_read_connect__name__[];
extern char pypylith3d_read_connect__doc__[];
extern "C"
PyObject * pypylith3d_read_connect(PyObject *, PyObject *);

// read coordinates
extern char pypylith3d_read_coords__name__[];
extern char pypylith3d_read_coords__doc__[];
extern "C"
PyObject * pypylith3d_read_coords(PyObject *, PyObject *);

// read differential forces
extern char pypylith3d_read_diff__name__[];
extern char pypylith3d_read_diff__doc__[];
extern "C"
PyObject * pypylith3d_read_diff(PyObject *, PyObject *);

// read time steps at which full output is desired
extern char pypylith3d_read_fuldat__name__[];
extern char pypylith3d_read_fuldat__doc__[];
extern "C"
PyObject * pypylith3d_read_fuldat(PyObject *, PyObject *);

// read time history info
extern char pypylith3d_read_hist__name__[];
extern char pypylith3d_read_hist__doc__[];
extern "C"
PyObject * pypylith3d_read_hist(PyObject *, PyObject *);

// read material history info
// these are not being used for now
// extern char pypylith3d_read_mathist__name__[];
// extern char pypylith3d_read_mathist__doc__[];
// extern "C"
// PyObject * pypylith3d_read_mathist(PyObject *, PyObject *);

// read prestresses
// extern char pypylith3d_read_prestr__name__[];
// extern char pypylith3d_read_prestr__doc__[];
// extern "C"
// PyObject * pypylith3d_read_prestr(PyObject *, PyObject *);

// read local coordinate rotations
extern char pypylith3d_read_skew__name__[];
extern char pypylith3d_read_skew__doc__[];
extern "C"
PyObject * pypylith3d_read_skew(PyObject *, PyObject *);

// read slippery node definitions
extern char pypylith3d_read_slip__name__[];
extern char pypylith3d_read_slip__doc__[];
extern "C"
PyObject * pypylith3d_read_slip(PyObject *, PyObject *);

// read split node definitions
extern char pypylith3d_read_split__name__[];
extern char pypylith3d_read_split__doc__[];
extern "C"
PyObject * pypylith3d_read_split(PyObject *, PyObject *);

// read state output information
extern char pypylith3d_read_stateout__name__[];
extern char pypylith3d_read_stateout__doc__[];
extern "C"
PyObject * pypylith3d_read_stateout(PyObject *, PyObject *);

// read time step group informations
extern char pypylith3d_read_timdat__name__[];
extern char pypylith3d_read_timdat__doc__[];
extern "C"
PyObject * pypylith3d_read_timdat(PyObject *, PyObject *);

// read traction boundary conditions
extern char pypylith3d_read_tractions__name__[];
extern char pypylith3d_read_tractions__doc__[];
extern "C"
PyObject * pypylith3d_read_tractions(PyObject *, PyObject *);

// read winkler force information
extern char pypylith3d_read_wink__name__[];
extern char pypylith3d_read_wink__doc__[];
extern "C"
PyObject * pypylith3d_read_wink(PyObject *, PyObject *);


// sort elements into element families
extern char pypylith3d_sort_elements__name__[];
extern char pypylith3d_sort_elements__doc__[];
extern "C"
PyObject * pypylith3d_sort_elements(PyObject *, PyObject *);

// Sort slippery nodes according to reordered elements
extern char pypylith3d_sort_slip_nodes__name__[];
extern char pypylith3d_sort_slip_nodes__doc__[];
extern "C"
PyObject * pypylith3d_sort_slip_nodes(PyObject *, PyObject *);

// Sort split nodes according to reordered elements
extern char pypylith3d_sort_split_nodes__name__[];
extern char pypylith3d_sort_split_nodes__doc__[];
extern "C"
PyObject * pypylith3d_sort_split_nodes(PyObject *, PyObject *);


// compute maximum number of nonzero entries in stiffness matrix
extern char pypylith3d_cmp_stiffsz__name__[];
extern char pypylith3d_cmp_stiffsz__doc__[];
extern "C"
PyObject * pypylith3d_cmp_stiffsz(PyObject *, PyObject *);

// create linked list for sparse matrix
extern char pypylith3d_lnklst__name__[];
extern char pypylith3d_lnklst__doc__[];
extern "C"
PyObject * pypylith3d_lnklst(PyObject *, PyObject *);

// localize id array for reference by element
extern char pypylith3d_local__name__[];
extern char pypylith3d_local__doc__[];
extern "C"
PyObject * pypylith3d_local(PyObject *, PyObject *);

// localize nfault array for reference by element
extern char pypylith3d_localf__name__[];
extern char pypylith3d_localf__doc__[];
extern "C"
PyObject * pypylith3d_localf(PyObject *, PyObject *);

// localize idx array for reference by element
extern char pypylith3d_localx__name__[];
extern char pypylith3d_localx__doc__[];
extern "C"
PyObject * pypylith3d_localx(PyObject *, PyObject *);

// create sparse matrix in modified sparse row format
extern char pypylith3d_makemsr__name__[];
extern char pypylith3d_makemsr__doc__[];
extern "C"
PyObject * pypylith3d_makemsr(PyObject *, PyObject *);


// drive time-dependent solution
extern char pypylith3d_viscos__name__[];
extern char pypylith3d_viscos__doc__[];
extern "C"
PyObject * pypylith3d_viscos(PyObject *, PyObject *);

// setup time-dependent solution
extern char pypylith3d_viscos_setup__name__[];
extern char pypylith3d_viscos_setup__doc__[];
extern "C"
PyObject * pypylith3d_viscos_setup(PyObject *, PyObject *);

// drive time-dependent solution
extern char pypylith3d_viscos_step__name__[];
extern char pypylith3d_viscos_step__doc__[];
extern "C"
PyObject * pypylith3d_viscos_step(PyObject *, PyObject *);

// cleanup time-dependent solution
extern char pypylith3d_viscos_cleanup__name__[];
extern char pypylith3d_viscos_cleanup__doc__[];
extern "C"
PyObject * pypylith3d_viscos_cleanup(PyObject *, PyObject *);


// write boundary conditions
extern char pypylith3d_write_bc__name__[];
extern char pypylith3d_write_bc__doc__[];
extern "C"
PyObject * pypylith3d_write_bc(PyObject *, PyObject *);

// write element node array
extern char pypylith3d_write_connect__name__[];
extern char pypylith3d_write_connect__doc__[];
extern "C"
PyObject * pypylith3d_write_connect(PyObject *, PyObject *);

// write nodal coordinates
extern char pypylith3d_write_coords__name__[];
extern char pypylith3d_write_coords__doc__[];
extern "C"
PyObject * pypylith3d_write_coords(PyObject *, PyObject *);

// write differential forces
extern char pypylith3d_write_diff__name__[];
extern char pypylith3d_write_diff__doc__[];
extern "C"
PyObject * pypylith3d_write_diff(PyObject *, PyObject *);

// write out element and prestress information
extern char pypylith3d_write_element_info__name__[];
extern char pypylith3d_write_element_info__doc__[];
extern "C"
PyObject * pypylith3d_write_element_info(PyObject *, PyObject *);

// write time steps at which full output is desired
extern char pypylith3d_write_fuldat__name__[];
extern char pypylith3d_write_fuldat__doc__[];
extern "C"
PyObject * pypylith3d_write_fuldat(PyObject *, PyObject *);

// write out global control information
extern char pypylith3d_write_global_info__name__[];
extern char pypylith3d_write_global_info__doc__[];
extern "C"
PyObject * pypylith3d_write_global_info(PyObject *, PyObject *);

// write load histories
extern char pypylith3d_write_hist__name__[];
extern char pypylith3d_write_hist__doc__[];
extern "C"
PyObject * pypylith3d_write_hist(PyObject *, PyObject *);

// write out material property info
extern char pypylith3d_write_props__name__[];
extern char pypylith3d_write_props__doc__[];
extern "C"
PyObject * pypylith3d_write_props(PyObject *, PyObject *);

// write local coordinate rotations
extern char pypylith3d_write_skew__name__[];
extern char pypylith3d_write_skew__doc__[];
extern "C"
PyObject * pypylith3d_write_skew(PyObject *, PyObject *);

// write slippery node entries
extern char pypylith3d_write_slip__name__[];
extern char pypylith3d_write_slip__doc__[];
extern "C"
PyObject * pypylith3d_write_slip(PyObject *, PyObject *);

// write out sparse matrix information
extern char pypylith3d_write_sparse_info__name__[];
extern char pypylith3d_write_sparse_info__doc__[];
extern "C"
PyObject * pypylith3d_write_sparse_info(PyObject *, PyObject *);

// write split node entries
extern char pypylith3d_write_split__name__[];
extern char pypylith3d_write_split__doc__[];
extern "C"
PyObject * pypylith3d_write_split(PyObject *, PyObject *);

// write split node entries for plot output
extern char pypylith3d_write_split_plot__name__[];
extern char pypylith3d_write_split_plot__doc__[];
extern "C"
PyObject * pypylith3d_write_split_plot(PyObject *, PyObject *);

// write state output information
extern char pypylith3d_write_stateout__name__[];
extern char pypylith3d_write_stateout__doc__[];
extern "C"
PyObject * pypylith3d_write_stateout(PyObject *, PyObject *);

// write out stress computation information
extern char pypylith3d_write_strscomp__name__[];
extern char pypylith3d_write_strscomp__doc__[];
extern "C"
PyObject * pypylith3d_write_strscomp(PyObject *, PyObject *);

// write out subiteration convergence information
extern char pypylith3d_write_subiter__name__[];
extern char pypylith3d_write_subiter__doc__[];
extern "C"
PyObject * pypylith3d_write_subiter(PyObject *, PyObject *);

// write time step data
extern char pypylith3d_write_timdat__name__[];
extern char pypylith3d_write_timdat__doc__[];
extern "C"
PyObject * pypylith3d_write_timdat(PyObject *, PyObject *);

// write traction BC
extern char pypylith3d_write_tractions__name__[];
extern char pypylith3d_write_tractions__doc__[];
extern "C"
PyObject * pypylith3d_write_tractions(PyObject *, PyObject *);

// write mesh info to UCD file
extern char pypylith3d_write_ucd_mesh__name__[];
extern char pypylith3d_write_ucd_mesh__doc__[];
extern "C"
PyObject * pypylith3d_write_ucd_mesh(PyObject *, PyObject *);

// write winkler BC
extern char pypylith3d_write_wink__name__[];
extern char pypylith3d_write_wink__doc__[];
extern "C"
PyObject * pypylith3d_write_wink(PyObject *, PyObject *);

// write slippery winkler BC
extern char pypylith3d_write_winkx__name__[];
extern char pypylith3d_write_winkx__doc__[];
extern "C"
PyObject * pypylith3d_write_winkx(PyObject *, PyObject *);


// scan boundary condition file
extern char pypylith3d_scan_bc__name__[];
extern char pypylith3d_scan_bc__doc__[];
extern "C"
PyObject * pypylith3d_scan_bc(PyObject *, PyObject *);

// scan connectivity file
extern char pypylith3d_scan_connect__name__[];
extern char pypylith3d_scan_connect__doc__[];
extern "C"
PyObject * pypylith3d_scan_connect(PyObject *, PyObject *);

// scan coordinates file
extern char pypylith3d_scan_coords__name__[];
extern char pypylith3d_scan_coords__doc__[];
extern "C"
PyObject * pypylith3d_scan_coords(PyObject *, PyObject *);

// scan differential forces file
extern char pypylith3d_scan_diff__name__[];
extern char pypylith3d_scan_diff__doc__[];
extern "C"
PyObject * pypylith3d_scan_diff(PyObject *, PyObject *);

// scan time step output info file
extern char pypylith3d_scan_fuldat__name__[];
extern char pypylith3d_scan_fuldat__doc__[];
extern "C"
PyObject * pypylith3d_scan_fuldat(PyObject *, PyObject *);

// scan time history definition file
extern char pypylith3d_scan_hist__name__[];
extern char pypylith3d_scan_hist__doc__[];
extern "C"
PyObject * pypylith3d_scan_hist(PyObject *, PyObject *);

// scan prestress file
// extern char pypylith3d_scan_prestr__name__[];
// extern char pypylith3d_scan_prestr__doc__[];
// extern "C"
// PyObject * pypylith3d_scan_prestr(PyObject *, PyObject *);

// scan local coordinate rotations file
extern char pypylith3d_scan_skew__name__[];
extern char pypylith3d_scan_skew__doc__[];
extern "C"
PyObject * pypylith3d_scan_skew(PyObject *, PyObject *);

// scan slippery node definitions file
extern char pypylith3d_scan_slip__name__[];
extern char pypylith3d_scan_slip__doc__[];
extern "C"
PyObject * pypylith3d_scan_slip(PyObject *, PyObject *);

// scan split node definitions file
extern char pypylith3d_scan_split__name__[];
extern char pypylith3d_scan_split__doc__[];
extern "C"
PyObject * pypylith3d_scan_split(PyObject *, PyObject *);

// scan time step group info file
extern char pypylith3d_scan_timdat__name__[];
extern char pypylith3d_scan_timdat__doc__[];
extern "C"
PyObject * pypylith3d_scan_timdat(PyObject *, PyObject *);

// scan traction boundary conditions file
extern char pypylith3d_scan_tractions__name__[];
extern char pypylith3d_scan_tractions__doc__[];
extern "C"
PyObject * pypylith3d_scan_tractions(PyObject *, PyObject *);

// scan winkler forces info file
extern char pypylith3d_scan_wink__name__[];
extern char pypylith3d_scan_wink__doc__[];
extern "C"
PyObject * pypylith3d_scan_wink(PyObject *, PyObject *);

// scan winkler forces info file for slippery nodes
extern char pypylith3d_scan_winkx__name__[];
extern char pypylith3d_scan_winkx__doc__[];
extern "C"
PyObject * pypylith3d_scan_winkx(PyObject *, PyObject *);

// Initialize material model info
extern char pypylith3d_matmod_def__name__[];
extern char pypylith3d_matmod_def__doc__[];
extern "C"
PyObject * pypylith3d_matmod_def(PyObject *, PyObject *);

// Precompute shape function info
extern char pypylith3d_preshape__name__[];
extern char pypylith3d_preshape__doc__[];
extern "C"
PyObject * pypylith3d_preshape(PyObject *, PyObject *);

// Precompute shape function info for element faces
extern char pypylith3d_preshape2d__name__[];
extern char pypylith3d_preshape2d__doc__[];
extern "C"
PyObject * pypylith3d_preshape2d(PyObject *, PyObject *);

// try_binio
extern char pypylith3d_try_binio__name__[];
extern char pypylith3d_try_binio__doc__[];
extern "C"
PyObject * pypylith3d_try_binio(PyObject *, PyObject *);


#endif // pypylith3d_libpylith3d_h

// end of file
