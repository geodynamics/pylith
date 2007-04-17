// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesive.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology::create()

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesive::FaultCohesive(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesive::~FaultCohesive(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::faults::FaultCohesive::FaultCohesive(const FaultCohesive& f) :
  Fault(f)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Adjust mesh topology for fault implementation.
void
pylith::faults::FaultCohesive::adjustTopology(ALE::Obj<ALE::Mesh>* mesh) const
{ // adjustTopology
  assert(0 != mesh);
  assert("" != label());

  // Get group of vertices associated with fault
  const ALE::Obj<int_section_type>& groupField = 
    (*mesh)->getIntSection(label());
  assert(!groupField.isNull());
  const int_section_type::chart_type& chart = groupField->getChart();

  // Create set with vertices on fault
  std::set<Mesh::point_type> points; // Vertices on fault
  const int numCells = (*mesh)->heightStratum(0)->size();
  for(int_section_type::chart_type::iterator c_iter = chart.begin();
      c_iter != chart.end();
      ++c_iter) {
    assert(!(*mesh)->depth(*c_iter));
    points.insert(*c_iter);
  } // for
  CohesiveTopology::create(*mesh, points);
} // adjustTopology


// End of file 
