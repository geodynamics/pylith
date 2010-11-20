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
#include "pylith/utils/array.hh" // HASA int_array

#include "IField.hh" // ISA ISection

// ISectionSpaces -------------------------------------------------------
/// Extend ALE::ISection to include spaces.
template<typename Point_, 
	 typename Value_, 
	 typename Alloc_ =ALE::malloc_allocator<Value_> >
class pylith::topology::ISectionSpaces : 
  public ALE::ISection<Point_, Value_, Alloc_>
{ // ISectionFields

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

  typedef ALE::ISection<Point_, Value_, Alloc_> base;
  typedef typename base::point_type point_type;
  typedef typename base::value_type value_type;
  typedef typename base::alloc_type alloc_type;
  typedef typename base::index_type index_type;
  typedef typename base::atlas_type atlas_type;
  typedef typename base::chart_type chart_type;

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

  /** Return only the values associated to this point, not its closure
   *
   * @param p Point associated with values.
   * @returns Values at point.
   */
  const value_type* restrictPoint(const point_type& p);

  /** Return only the values associated to this point, not its closure
   *
   * @param p Point associated with values.
   * @param values Array in which to store values.
   * @param size Size of array.
   */
  void restrictPoint(const point_type& p,
		     value_type* const values,
		     const int size);

  /** Get number of spaces.
   *
   * @returns Number of spaces.
   */
  int getNumSpaces(void) const;

  /// Add space.
  void addSpace();
  
  /** Set fiber dimension of space.
   *
   * @param space Space identifier (index).
   * @param fiberDim Fiberdimension of space.
   */
  void spaceFiberDimension(const int space,
			   const int fiberDim);
  
  /** Get field associated with space.
   *
   * @param fibration Index of space.
   */
  ALE::Obj<ALE::IGeneralSection<Point_, Value_, Alloc_> >
    getFibration(const int space);
  
// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  int_array _spaces; ///< Fiber dimensions of spaces.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  ISectionSpaces(const ISectionSpaces&); ///< Not implemented
  const ISectionSpaces& operator=(const ISectionSpaces&); ///< Not implemented

}; // ISectionSpaces

#include "ISectionSpaces.cc"

#endif // pylith_topology_isectionspaces_hh


// End of file 
