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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/problems/Implicit.hh
 *
 * @brief Object for implicit time integration.
 */

#if !defined(pylith_problems_implicit_hh)
#define pylith_problems_implicit_hh

// Include directives ---------------------------------------------------
#include "Formulation.hh" // ISA Formulation

// Implicit ---------------------------------------------------------
/** @brief Object for implicit time integration.
 *
 * Implicit time stepping associated with quasi-static problems.
 */

class pylith::problems::Implicit : public Formulation
{ // Implicit
  friend class TestImplicit; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Implicit(void);

  /// Destructor
  ~Implicit(void);

  /// Compute rate fields (velocity and/or acceleration) at time t.
  void calcRateFields(void);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Implicit(const Implicit&); ///< Not implemented
  const Implicit& operator=(const Implicit&); ///< Not implemented

}; // Implicit

#endif // pylith_problems_implicit_hh


// End of file 
