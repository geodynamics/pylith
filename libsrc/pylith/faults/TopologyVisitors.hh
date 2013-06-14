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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/TopologyVisitors.hh
 *
 * @brief C++ objects for visitors to topology data structures.
 */

#if !defined(pylith_faults_topologyvisitors_hh)
#define pylith_faults_topologyvisitors_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

// ReplaceVisitor -------------------------------------------------------
/// Visitor for replacement.
template<typename Sieve, typename Renumbering>
class pylith::faults::ReplaceVisitor {
private:
  typedef typename Sieve::point_type point_type;
protected:
  Renumbering& renumbering;
  const int    size;
  int          i;
  const int    debug;
  point_type  *points;
  bool         mapped;
public:
  ReplaceVisitor(Renumbering& r,
		 const int size,
		 const int debug =0);
  ~ReplaceVisitor(void);
  void visitPoint(const point_type& point);
  void visitArrow(const typename Sieve::arrow_type&);
  const point_type *getPoints(void);
  bool mappedPoint(void);
  void clear(void);
};

// ClassifyVisitor ------------------------------------------------------
/// Visitor for classification.
template<typename Sieve>
class pylith::faults::ClassifyVisitor {
public:
  typedef typename Sieve::point_type point_type;
  typedef std::set<point_type> PointSet;
protected:
  const Sieve&     sieve;
  const PointSet&  replaceCells;
  const PointSet&  noReplaceCells;
  const point_type firstCohesiveCell;
  const int        faceSize;
  const int        debug;
  PointSet         vReplaceCells;
  PointSet         vNoReplaceCells;
  bool             modified;
  bool             setupMode;
  int              size;
  ALE::ISieveVisitor::PointRetriever<Sieve> pR;
public:
  ClassifyVisitor(const Sieve& s,
		  const PointSet& rC,
		  const PointSet& nrC,
		  const point_type& fC,
		  const int fS,
		  const int debug =0);
  ~ClassifyVisitor(void);
  void visitPoint(const point_type& point);
  void visitArrow(const typename Sieve::arrow_type&);
  const PointSet& getReplaceCells() const;
  const PointSet& getNoReplaceCells() const;
  bool getModified() const;
  int getSize() const;
  void setMode(const bool isSetup);
  void reset();
};

#include "TopologyVisitors.cc" // template definitions

#endif // pylith_faults_topologyvisitors_hh


// End of file 
