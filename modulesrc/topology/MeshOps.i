// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================

/**
 * @file modulesrc/topology/MeshOps.i
 *
 * @brief Python interface to C++ MeshOps.
 */

%inline %{
  /** Nondimensionalize the finite-element mesh.
   *
   * @param mesh Finite-element mesh.
   * @param normalizer Nondimensionalizer.
   */
  void
  MeshOps_nondimensionalize(pylith::topology::Mesh* const mesh,
			    const spatialdata::units::Nondimensional& normalizer) {
    pylith::topology::MeshOps::nondimensionalize(mesh, normalizer);
  } // nondimensionalize
%}

// End of file
