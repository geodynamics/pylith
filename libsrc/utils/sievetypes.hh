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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/utils/sievetypes.hh
 *
 * @brief Aliases for Sieve types.
 */

#if !defined(pylith_utils_sievetypes_hh)
#define pylith_utils_sievetypes_hh

#include <petscdmmesh.hh> // PETSc Mesh

namespace pylith {

  /// Sieve mesh (default, fast access with set sizes).
  typedef ALE::IMesh<PetscInt,PetscScalar> SieveMesh;

  /// Sieve mesh (flexible, slower access without set sizes).
  typedef ALE::Mesh<PetscInt,PetscScalar> SieveFlexMesh;

  /// Sieve submesh.
  typedef ALE::IMesh<PetscInt,PetscScalar,ALE::LabelSifter<int, SieveMesh::point_type> > SieveSubMesh;

} // pylith

#endif // pylith_utils_sievetypes_hh


// End of file
