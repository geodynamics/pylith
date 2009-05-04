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

/** @file libsrc/bc/bcfwd.hh
 *
 * @brief Forward declarations for PyLith boundary condition objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_bc_bcfwd_hh)
#define pylith_bc_bcfwd_hh

namespace pylith {
  namespace bc {

    class BoundaryCondition;
    class DirichletBC;
    class DirichletBoundary;
    class Neumann;
    class AbsorbingDampers;

  } // bc
} // pylith


#endif // pylith_bc_bcfwd_hh


// End of file 
