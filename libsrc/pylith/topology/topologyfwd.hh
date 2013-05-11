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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

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

namespace ALE {
  template<typename cellrefiner_type> class MeshRefiner;
  class RefineEdges2;
  class CellRefinerTri3;
  class CellRefinerTet4;

  class RefineFace4Edges2;
  class CellRefinerQuad4;

  class RefineVol8Face4Edges2;
  class CellRefinerHex8;

  class MeshOrder;
} // ALE

namespace pylith {
  namespace topology {

    class Mesh;
    class SubMesh;
    class MeshOps;
    class CoordsVisitor;
    class SubMeshIS;
    class Stratum;
    class StratumIS;
    
    class FieldBase;
    template<typename mesh_type> class Field;
    template<typename field_type> class Fields;
    class VecVisitorMesh;
    class VecVisitorSubMesh;

    class SolutionFields;

    class Jacobian;
    class MatVisitorMesh;
    class MatVisitorSubMesh;

    class Distributor;

    class RefineUniform;

    class ReverseCuthillMcKee;

  } // topology
} // pylith


#endif // pylith_topology_topologyfwd_hh


// End of file 
