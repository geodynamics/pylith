// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/topology/ISectionSpaces.hh
 *
 * @brief Extend ALE::ISection to include spaces.
 */

#if !defined(pylith_topology_isectionspaces_hh)
#define pylith_topology_isectionspaces_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "IField.hh"

// ISectionSpaces -------------------------------------------------------
/// Extend ALE::ISection to include spaces.
template<typename Point_, 
	 typename Value_, 
	 typename Alloc_ =ALE::malloc_allocator<Value_> >
class pylith::topology::ISectionSpaces : 
  public ALE::ISection<Point_, Value_, Alloc_> {
{ // ISectionFields

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param comm MPI communicator.
   * @param debug Debugging flag.
   */
  ISectionSpaces(MPI_Comm comm,
		 const int debug =0);

  /** Constructor with communicator and max/min for chart.
   *
   * @param comm MPI communicator.
   * @param min Minimum point for chart.
   * @param max Maximum point for chart.
   * @param debug Debugging flag.
   */ 
  ISectionSpaces(MPI_Comm comm,
		 const point_type& min,
		 const point_type& max,
		 const int debug =0);

  /** Constructor with atlas.
   *
   * @param atlas Atlas for section.
   */
  ISectionSpaces(const ALE::Obj<atlas_type>& atlas);

  /// Destructor.
  virtual
  ~ISectionSpaces(void);

  /** Get number of spaces.
   *
   * @returns Number of spaces.
   */
  int getNumSpaces(void) const;

  /// Add space.
  void addSpace(void);
  
  /** Get field associated with space.
   *
   * @param fibration Index of space.
   */
  ALE::Obj<IGeneralSection<point_type, value_type, alloc_type> >&
    getFibration(const int space);
  
// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  map_type _fields; ///< Fields without constraints over a common set of points.
  ALE::Obj<section_type> _section; ///< Section containing fields.
  const mesh_type& _mesh; ///< Mesh associated with fields.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  ISectionSpaces(const ISectionSpaces&); ///< Not implemented
  const ISectionSpaces& operator=(const ISectionSpaces&); ///< Not implemented

}; // ISectionSpaces

#include "ISectionSpaces.cc"

#endif // pylith_topology_isectionspaces_hh


// End of file 
