// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::ISectionSpaces(MPI_Comm comm,
	       const int debug =0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Constructor with communicator and max/min for chart.
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::ISectionSpaces(MPI_Comm comm,
		 const point_type& min,
		 const point_type& max,
		 const int debug =0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Constructor with atlas.
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::ISectionSpaces(const ALE::Obj<atlas_type>& atlas)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::~ISectionSpaces(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get number of spaces.
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
int 
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::getNumSpaces(void) const
{ // getNumSpaces
} // getNumSpaces

// ----------------------------------------------------------------------
// Add space.
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void 
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::addSpace(void)
{ // addSpace
} // addSPace
  
// ----------------------------------------------------------------------
// Get field associated with space.
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
ALE::Obj<IGeneralSection<point_type, value_type, alloc_type> >&
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::getFibration(const int space)
{ // getFibration
} // getFibration


// End of file 
