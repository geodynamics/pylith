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

/** @file libsrc/problems/problemsfwd.hh
 *
 * @brief Forward declarations for PyLith problems objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_problems_problemsfwd_hh)
#define pylith_problems_problemsfwd_hh

namespace pylith {
  namespace problems {

    class Formulation;

    class Solver;
    class SolverLinear;
    class SolverNonlinear;

  } // problems
} // pylith


#endif // pylith_problems_problemsfwd_hh


// End of file 
