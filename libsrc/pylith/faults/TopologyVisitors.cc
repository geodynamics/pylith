// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

// ----------------------------------------------------------------------
template<typename Sieve, typename Renumbering>
pylith::faults::ReplaceVisitor<Sieve,Renumbering>::ReplaceVisitor(
							  Renumbering& r,
							  const int size,
							  const int debug) :
  renumbering(r),
  size(size),
  i(0),
  debug(debug)
{ // constructor
  this->points = new point_type[this->size];
  this->mapped = false;
} // constructor

// ----------------------------------------------------------------------
template<typename Sieve, typename Renumbering>
pylith::faults::ReplaceVisitor<Sieve,Renumbering>::~ReplaceVisitor()
{ // destructor
  delete[] this->points;
} // destructor

// ----------------------------------------------------------------------
template<typename Sieve, typename Renumbering>
void
pylith::faults::ReplaceVisitor<Sieve,Renumbering>::visitPoint(
						     const point_type& point)
{ // visitPoint
  if (i >= this->size)
    throw ALE::Exception("Too many points for ReplaceVisitor");
  if (this->renumbering.find(point) != this->renumbering.end()) {
    if (debug)
      std::cout << "    point " << this->renumbering[point] << std::endl;
    points[i] = this->renumbering[point];
    this->mapped = true;
  } else {
    if (debug) std::cout << "    point " << point << std::endl;
    points[i] = point;
  } // if/else
  ++i;
} // visitPoint

// ----------------------------------------------------------------------
template<typename Sieve, typename Renumbering>
void
pylith::faults::ReplaceVisitor<Sieve,Renumbering>::visitArrow(
					   const typename Sieve::arrow_type&)
{ // visitArrow
} // visitArrow

// ----------------------------------------------------------------------
template<typename Sieve, typename Renumbering>
inline
const typename Sieve::point_type*
pylith::faults::ReplaceVisitor<Sieve,Renumbering>::getPoints(void)
{ // getPoints
  return this->points;
} // getPoints

// ----------------------------------------------------------------------
template<typename Sieve, typename Renumbering>
inline
bool
pylith::faults::ReplaceVisitor<Sieve,Renumbering>::mappedPoint(void)
{ // mappedPoint
  return this->mapped;
} // mappedPoint

// ----------------------------------------------------------------------
template<typename Sieve, typename Renumbering>
inline
void
pylith::faults::ReplaceVisitor<Sieve,Renumbering>::clear(void)
{ // clear
  this->i = 0; this->mapped = false;
} // clear


// ClassifyVisitor ------------------------------------------------------
template<typename Sieve>
pylith::faults::ClassifyVisitor<Sieve>::ClassifyVisitor(const Sieve& s,
							const PointSet& rC,
							const PointSet& nrC,
							const point_type& fC,
							const int fS,
							const int debug) :
  sieve(s),
  replaceCells(rC),
  noReplaceCells(nrC),
  firstCohesiveCell(fC),
  faceSize(fS),
  debug(debug),
  modified(false),
  setupMode(true),
  size(0)
{ // constructor
  pR.setSize(s.getMaxConeSize());
} // constructor

// ----------------------------------------------------------------------
template<typename Sieve>
pylith::faults::ClassifyVisitor<Sieve>::~ClassifyVisitor(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
template<typename Sieve>
void
pylith::faults::ClassifyVisitor<Sieve>::visitPoint(const point_type& point)
{ // visitPoint
  if (this->setupMode) {
    if (replaceCells.find(point) != replaceCells.end())
      vReplaceCells.insert(point);
    if (noReplaceCells.find(point) != noReplaceCells.end())
      vNoReplaceCells.insert(point);
    if (point >= firstCohesiveCell) return;
    this->modified = true;
    this->size++;
    return;
  } // if
  bool classified = false;
    
  if (debug) {
    std::cout << "Checking neighbor " << point << std::endl;
    ALE::ISieveVisitor::PointRetriever<Sieve> pR(9);

    sieve.cone(point, pR);
    for(int i = 0; i < pR.getSize(); ++i) {
      std::cout << "  cone point " << pR.getPoints()[i] << std::endl;
    }
  }
  if (vReplaceCells.find(point) != vReplaceCells.end()) {
    if (debug) 
      std::cout << "  already in replaceCells" << std::endl;
    return;
  } // if
  if (vNoReplaceCells.find(point) != vNoReplaceCells.end()) {
    if (debug) 
      std::cout << "  already in noReplaceCells" << std::endl;
    return;
  } // if
  if (point >= firstCohesiveCell) {
    if (debug) 
      std::cout << "  already a cohesive cell" << std::endl;
    return;
  } // if
    // If neighbor shares a face with anyone in replaceCells, then add
  for (typename PointSet::const_iterator c_iter = vReplaceCells.begin();
      c_iter != vReplaceCells.end();
      ++c_iter) {
    sieve.meet(*c_iter, point, pR);
      
    if (pR.getSize() == faceSize) {
      if (debug)
	std::cout << "    Scheduling " << point << " for replacement"
		  << std::endl;
      vReplaceCells.insert(point);
      modified   = true;
      classified = true;
      pR.clear();
      break;
    } // if
    pR.clear();
  } // for
  if (classified)
    return;
  // It is unclear whether taking out the noReplace cells will speed this up
  for (typename PointSet::const_iterator c_iter = vNoReplaceCells.begin();
      c_iter != vNoReplaceCells.end();
      ++c_iter) {
    sieve.meet(*c_iter, point, pR);
      
    if (pR.getSize() == faceSize) {
      if (debug) 
	std::cout << "    Scheduling " << point << " for no replacement"
		  << std::endl;
      vNoReplaceCells.insert(point);
      modified   = true;
      classified = true;
      pR.clear();
      break;
    } // for
    pR.clear();
  } // for
} // visitPoint

// ----------------------------------------------------------------------
template<typename Sieve>
inline
void
pylith::faults::ClassifyVisitor<Sieve>::visitArrow(const typename Sieve::arrow_type&)
{ // visitArrow
} // visitArrow

// ----------------------------------------------------------------------
template<typename Sieve>
inline
const std::set<typename Sieve::point_type>&
pylith::faults::ClassifyVisitor<Sieve>::getReplaceCells(void) const
{ // getReplaceCells
  return this->vReplaceCells;
} // getReplaceCells

// ----------------------------------------------------------------------
template<typename Sieve>
inline
const std::set<typename Sieve::point_type>&
pylith::faults::ClassifyVisitor<Sieve>::getNoReplaceCells() const
{ // getNoReplaceCells
  return this->vNoReplaceCells;
} // getNoReplaceCells

// ----------------------------------------------------------------------
template<typename Sieve>
inline
bool
pylith::faults::ClassifyVisitor<Sieve>::getModified() const
{ // getModified
  return this->modified;
} // getModified

// ----------------------------------------------------------------------
template<typename Sieve>
inline
int
pylith::faults::ClassifyVisitor<Sieve>::getSize() const
{ // getSize
  return this->size;
} // getSize

// ----------------------------------------------------------------------
template<typename Sieve>
inline
void
pylith::faults::ClassifyVisitor<Sieve>::setMode(const bool isSetup)
{ // setMode
  this->setupMode = isSetup;
} // setMode

// ----------------------------------------------------------------------
template<typename Sieve>
inline
void
pylith::faults::ClassifyVisitor<Sieve>::reset(void)
{ // reset
  this->modified = false;
} // reset


// End of file 
