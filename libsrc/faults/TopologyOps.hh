// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/TopologyOps.hh
 *
 * @brief C++ object to manage creation of cohesive cells.
 */

#if !defined(pylith_faults_topologyops_hh)
#define pylith_faults_topologyops_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

// TopologyOps ----------------------------------------------------------
class pylith::faults::TopologyOps
{ // class TopologyOps
public :
  typedef std::set<SieveMesh::point_type> PointSet;
  typedef std::vector<sieve_type::point_type> PointArray;
  typedef std::pair<sieve_type::point_type, int> oPoint_type;
  typedef std::vector<oPoint_type>  oPointArray;

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  template<class InputPoints>
  static
  bool compatibleOrientation(const topology::& mesh,
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
  void computeCensoredDepth(const ALE::Obj<Mesh::label_type>& depth,
			    const ALE::Obj<Mesh::sieve_type>& sieve,
			    const Mesh::point_type& firstCohesiveCell);
  
  static
  void classifyCells(const ALE::Obj<Mesh::sieve_type>& sieve,
		     const Mesh::point_type& vertex,
		     const int depth,
		     const int faceSize,
		     const Mesh::point_type& firstCohesiveCell,
		     PointSet& replaceCells,
		     PointSet& noReplaceCells,
		     const int debug);
  
  static
  void createFaultSieveFromVertices(const int dim,
				    const int firstCell,
				    const PointSet& faultVertices,
				    const ALE::Obj<Mesh>& mesh,
				    const ALE::Obj<ALE::Mesh::arrow_section_type>& orientation,
				    const ALE::Obj<ALE::Mesh::sieve_type>& faultSieve,
				    const bool flipFault);
  
  static
  void createFaultSieveFromFaces(const int dim,
				 const int firstCell,
				 const int numFaces,
				 const int faultVertices[],
				 const int faultCells[],
				 const ALE::Obj<Mesh>& mesh,
				 const ALE::Obj<ALE::Mesh::arrow_section_type>& orientation,
				 const ALE::Obj<ALE::Mesh::sieve_type>& faultSieve);

  static
  void orientFaultSieve(const int dim,
			const ALE::Obj<Mesh>& mesh,
			const ALE::Obj<ALE::Mesh::arrow_section_type>& orientation,
			const ALE::Obj<ALE::Mesh>& fault);
}; // class CohesiveTopology

#endif // pylith_faults_cohesivetopology_hh


// End of file 
