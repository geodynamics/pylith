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
 * @file libsrc/problems/Explicit.hh
 *
 * @brief Object for explicit time integration.
 */

#if !defined(pylith_problems_explicit_hh)
#define pylith_problems_explicit_hh

// Include directives ---------------------------------------------------
#include "Formulation.hh" // ISA Formulation

// Explicit ---------------------------------------------------------
/** @brief Object for explicit time integration.
 *
 * Explicit time stepping associated with dynamic problems.
 */

class pylith::problems::Explicit : public Formulation
{ // Explicit
  friend class TestExplicit; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Explicit(void);

  /// Destructor
  ~Explicit(void);

  /// Compute rate fields (velocity and/or acceleration) at time t.
  void calcRateFields(void);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Explicit(const Explicit&); ///< Not implemented
  const Explicit& operator=(const Explicit&); ///< Not implemented

}; // Explicit

#endif // pylith_problems_explicit_hh


// End of file 
