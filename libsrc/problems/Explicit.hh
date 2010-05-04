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

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  /// Setup rate fields.
  void _setupRateFields(void);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Explicit(const Explicit&); ///< Not implemented
  const Explicit& operator=(const Explicit&); ///< Not implemented

}; // Explicit

#endif // pylith_problems_explicit_hh


// End of file 
