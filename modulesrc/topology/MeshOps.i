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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

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
