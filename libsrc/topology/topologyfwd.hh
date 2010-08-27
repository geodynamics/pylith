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

/** @file libsrc/topology/topologyfwd.hh
 *
 * @brief Forward declarations for PyLith topology objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_topology_topologyfwd_hh)
#define pylith_topology_topologyfwd_hh

namespace pylith {
  namespace topology {

    class Mesh;
    class SubMesh;
    class MeshOps;

    class FieldBase;
    template<typename mesh_type> class Field;
    template<typename field_type> class Fields;
    template<typename field_type, int fiberDimTotal> class FieldsNew;
    class SolutionFields;

    class Jacobian;

    class Distributor;

    class MeshRefiner;
    class RefineUniform;

    class ReverseCuthillMcKee;

  } // topology
} // pylith


#endif // pylith_topology_topologyfwd_hh


// End of file 
