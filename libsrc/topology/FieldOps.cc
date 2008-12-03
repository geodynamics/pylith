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

#include <portinfo>

#include "FieldOps.hh" // implementation of class methods

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Copy values from one section to another.
void
pylith::topology::FieldOps::copyValues(const ALE::Obj<real_section_type>& dest,
				       const ALE::Obj<real_section_type>& src)
{ // copyValues
  typedef real_section_type::chart_type chart_type;

  assert(!dest.isNull());
  assert(!src.isNull());

  const chart_type& chartSrc = src->getChart();
  const chart_type& chartDest = dest->getChart();
  const chart_type::const_iterator chartEnd = chartSrc.end();
  for (chart_type::const_iterator c_iter = chartSrc.begin();
       c_iter != chartEnd;
       ++c_iter) {
    assert(dest->getFiberDimension(*c_iter) == 
	   src->getFiberDimension(*c_iter));
    dest->updatePoint(*c_iter, src->restrictPoint(*c_iter));
  } // for
} // copyValues


// End of file 
