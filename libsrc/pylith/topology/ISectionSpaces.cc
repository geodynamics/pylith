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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::ISectionSpaces(MPI_Comm comm,
										     const int debug) :
  ALE::ISection<point_type, value_type, alloc_type>(comm, debug)
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
										     const int debug) :
  ALE::ISection<point_type, value_type, alloc_type>(comm, min, max, debug)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Constructor with atlas.
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::ISectionSpaces(const ALE::Obj<atlas_type>& atlas) :
  ALE::ISection<point_type, value_type, alloc_type>(atlas)
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
// Return only the values associated to this point, not its closure
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
const typename pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::value_type*
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::restrictPoint(const point_type& p)
{ // restrictPoint
  return ALE::ISection<point_type, value_type, alloc_type>::restrictPoint(p);
} // restrictPoint

// ----------------------------------------------------------------------
// Return only the values associated to this point, not its closure
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::restrictPoint(const point_type& p,
		     value_type* const values,
		     const int size)
{ // restrictPoint
  assert(values);

  assert(this->_array);
  if (this->hasPoint(p)) {
    const int fiberDim = this->getFiberDimension(p);
    assert(size == fiberDim);
    const value_type* valuesPoint = this->restrictPoint(p);
    for (int i=0; i < fiberDim; ++i)
      values[i] = valuesPoint[i];
  } // if/else
} // restrictPoint

// ----------------------------------------------------------------------
// Get number of spaces.
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
int 
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::getNumSpaces(void) const
{ // getNumSpaces
  return _spaces.size();
} // getNumSpaces

// ----------------------------------------------------------------------
// Add space.
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void 
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::addSpace(void)
{ // addSpace
  _spaces.resize(getNumSpaces()+1);
} // addSPace
  
// ----------------------------------------------------------------------
// Set fiber dimension of space.
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::spaceFiberDimension(const int space,
											const int fiberDim)
{ // setFiberDimension
  assert(space >= 0 && space < _spaces.size());
  assert(fiberDim > 0);
  _spaces[space] = fiberDim;
} // setFiberDimension
  
// ----------------------------------------------------------------------
// Get field associated with space.
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
ALE::Obj<ALE::IGeneralSection<point_type, value_type, alloc_type> >
pylith::topology::ISectionSpaces<point_type, value_type, alloc_type>::getFibration(const int space)
{ // getFibration
  assert(space >= 0 && space < _spaces.size());

  typedef typename ALE::IGeneralSection<point_type, value_type, alloc_type>::chart_type IGeneralSection_chart_type;
  typedef typename ALE::IGeneralSection<point_type, value_type, alloc_type>::atlas_type IGeneralSection_atlas_type;

  ALE::Obj<ALE::IGeneralSection<point_type, value_type, alloc_type> > field =
    new ALE::IGeneralSection<point_type, value_type, alloc_type>(this->comm(),
							    this->debug());
  assert(!field.isNull());
  field->setChart(this->getChart());
  const chart_type& chart = this->getChart();

  // Copy sizes
  const int fiberDim = _spaces[space];
  const typename chart_type::const_iterator chartEnd = chart.end();
  for(typename chart_type::const_iterator c_iter = chart.begin();
      c_iter != chartEnd;
      ++c_iter) {
    if (this->getFiberDimension(*c_iter) > 0)
      field->setFiberDimension(*c_iter, fiberDim);
  } // for
  field->allocateStorage();
  ALE::Obj<IGeneralSection_atlas_type> newAtlas = 
    new IGeneralSection_atlas_type(this->comm(), this->debug());
  assert(!newAtlas.isNull());
  const IGeneralSection_chart_type& newChart = field->getChart();

  // Compute relative offset for space
  int offsetSpace = 0;
  for (int i=0; i < space; ++i)
    offsetSpace += _spaces[i];

  // Copy offsets
  const ALE::Obj<atlas_type>& thisAtlas = this->getAtlas();
  assert(!thisAtlas.isNull());
  newAtlas->setChart(newChart);
  newAtlas->allocatePoint();
  const typename IGeneralSection_chart_type::const_iterator newChartEnd =
    newChart.end();
  typename ALE::IGeneralSection<point_type, value_type, alloc_type>::index_type idx;
  for (typename IGeneralSection_chart_type::const_iterator c_iter=newChart.begin();
       c_iter != newChartEnd;
       ++c_iter) {
    if (this->getFiberDimension(*c_iter) > 0) {
      idx.prefix = fiberDim;
      const int c_offset = thisAtlas->restrictPoint(*c_iter)[0].index;
      idx.index  = c_offset + offsetSpace;
    } else {
      idx.prefix = 0;
      idx.index = 0;
    } // if/else
    newAtlas->addPoint(*c_iter);
    newAtlas->updatePoint(*c_iter, &idx);
  } // for
  field->replaceStorage(this->_array, true, this->sizeWithBC());
  field->setAtlas(newAtlas);

  return field;
} // getFibration


// End of file 
