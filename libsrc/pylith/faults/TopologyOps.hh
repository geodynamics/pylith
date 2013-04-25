// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/TopologyOps.hh
 *
 * @brief C++ helper object for creation of cohesive cells.
 */

#if !defined(pylith_faults_topologyops_hh)
#define pylith_faults_topologyops_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/sievetypes.hh" // USE SieveFlexMesh

// TopologyOps ----------------------------------------------------------
/// Helper object for creation of cohesive cells.
class pylith::faults::TopologyOps
{ // class TopologyOps

  // PUBLIC TYPEDEFS ////////////////////////////////////////////////////
public :
  typedef pylith::topology::Mesh::SieveMesh SieveMesh;
  typedef std::set<SieveMesh::point_type> PointSet;
  typedef std::vector<SieveMesh::sieve_type::point_type> PointArray;
  typedef std::pair<SieveMesh::sieve_type::point_type, int> oPoint_type;
  typedef std::vector<oPoint_type>  oPointArray;

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  template<class InputPoints>
  static
  bool compatibleOrientation(const ALE::Obj<SieveMesh>& mesh,
			     const SieveMesh::point_type& p,
			     const SieveMesh::point_type& q,
			     const int numFaultCorners,
			     const int faultFaceSize,
			     const int faultDepth,
			     const ALE::Obj<InputPoints>& points,
			     int indices[],
			     PointArray *origVertices,
			     PointArray *faceVertices,
			     PointArray *neighborVertices);

  static
  void computeCensoredDepth(const ALE::Obj<SieveMesh::label_type>& depth,
			    const ALE::Obj<SieveMesh::sieve_type>& sieve,
			    const SieveMesh::point_type& firstCohesiveCell);
  
  static
  void classifyCells(const ALE::Obj<SieveMesh::sieve_type>& sieve,
		     const SieveMesh::point_type& vertex,
		     const int depth,
		     const int faceSize,
		     const SieveMesh::point_type& firstCohesiveCell,
		     PointSet& replaceCells,
		     PointSet& noReplaceCells,
		     const int debug);
  static
  void classifyCellsDM(DM dmMesh,
		     PetscInt vertex,
		     const int depth,
		     const int faceSize,
		     PetscInt firstCohesiveCell,
		     PointSet& replaceCells,
		     PointSet& noReplaceCells,
		     const int debug);
  
  static
  void createFaultSieveFromVertices(const int dim,
				    const int firstCell,
				    const PointSet& faultVertices,
				    const ALE::Obj<SieveMesh>& mesh,
                    const ALE::Obj<SieveFlexMesh::arrow_section_type>& orientation,
				    const ALE::Obj<SieveFlexMesh::sieve_type>& faultSieve,
				    const bool flipFault);
  
  static
  void createFaultSieveFromFaces(const int dim,
				 const int firstCell,
				 const int numFaces,
				 const int faultVertices[],
				 const int faultCells[],
				 const ALE::Obj<SieveMesh>& mesh,
				 const ALE::Obj<SieveFlexMesh::arrow_section_type>& orientation,
				 const ALE::Obj<SieveFlexMesh::sieve_type>& faultSieve);

  static
  void orientFaultSieve(const int dim,
			const ALE::Obj<SieveMesh>& mesh,
			const ALE::Obj<SieveFlexMesh::arrow_section_type>& orientation,
			const ALE::Obj<SieveFlexMesh>& fault);
}; // class CohesiveTopology

#endif // pylith_faults_cohesivetopology_hh


// End of file 
