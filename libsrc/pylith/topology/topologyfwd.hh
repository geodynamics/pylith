// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file libsrc/topology/topologyfwd.hh
 *
 * @brief Forward declarations for PyLith topology objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_topology_topologyfwd_hh)
#define pylith_topology_topologyfwd_hh

#include "pylith/utils/types.hh"

namespace pylith {
    namespace topology {
        class Mesh;
        class MeshOps;
        class CoordsVisitor;
        class SubmeshIS;
        class Stratum;
        class StratumIS;

        class FieldBase;
        class Field;
        class VecVisitorMesh;
        class VecVisitorSubmesh;

        class FE;
        class FieldFactory;
        class FieldOps;
        class FieldQuery;

        class MatVisitorMesh;
        class MatVisitorSubmesh;

        class Distributor;
        class RefineUniform;
        class ReverseCuthillMcKee;

    } // topology
} // pylith

#endif // pylith_topology_topologyfwd_hh

// End of file
