// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


%inline %{
  /** Nondimensionalize the finite-element mesh.
   *
   * @param mesh Finite-element mesh.
   * @param scales Nondimensionalizer.
   */
  void
  MeshOps_nondimensionalize(pylith::topology::Mesh* const mesh,
			    const pylith::scales::Scales& scales) {
    pylith::topology::MeshOps::nondimensionalize(mesh, scales);
  } // nondimensionalize
%}


// End of file
