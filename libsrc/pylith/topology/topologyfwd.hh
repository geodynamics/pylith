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
// Copyright (c) 2010-2011 University of California, Davis
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

#include "pylith/utils/sievetypes.hh"
#include "pylith/utils/types.hh"

namespace ALE {
  template<typename point_type, 
	   typename value_type,
	   typename allocator> class IUniformSectionDS;

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

    class FieldBase;
    template<typename mesh_type, 
	     typename section_type =ALE::IGeneralSection<pylith::SieveMesh::point_type, PylithScalar> > class Field;
    template<typename field_type> class Fields;
    template<typename mesh_type> class FieldsNew;

    class SolutionFields;

    class Jacobian;

    class Distributor;

    class RefineUniform;

    class ReverseCuthillMcKee;

  template<typename point_type, 
	   typename value_type,
	   typename allocator> class ISectionSpaces;
  } // topology
} // pylith


#endif // pylith_topology_topologyfwd_hh


// End of file 
