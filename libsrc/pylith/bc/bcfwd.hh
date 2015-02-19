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
    class BoundaryConditionPoints;
    class BCIntegratorSubMesh;
    class TimeDependent;
    class TimeDependentPoints;
    class TimeDependentSubMesh;
    class DirichletBC;
    class DirichletBoundary;
    class Neumann;
    class Neumann_NEW;
    class AbsorbingDampers;
    class PointForce;

  } // bc
} // pylith


#endif // pylith_bc_bcfwd_hh


// End of file 
