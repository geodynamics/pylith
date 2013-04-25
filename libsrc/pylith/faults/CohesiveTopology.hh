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

/** @file libsrc/faults/CohesiveTopology.hh
 *
 * @brief C++ object to manage creation of cohesive cells.
 */

#if !defined(pylith_faults_cohesivetopology_hh)
#define pylith_faults_cohesivetopology_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/topology/Mesh.hh" // USES Mesh::IntSection
#include "pylith/utils/sievetypes.hh" // USE ALE::Obj and SieveFlexMesh

// CohesiveTopology -----------------------------------------------------
/// Creation of cohesive cells.
class pylith::faults::CohesiveTopology
{ // class CohesiveTopology

private :
  typedef pylith::topology::Mesh::SieveMesh::point_type point_type;

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Create the fault mesh.
   *
   * @param faultMesh Finite-element mesh of fault (output).
   * @param faultBoundary Finite-element mesh of fault boundary (output).
   * @param mesh Finite-element mesh of domain.
   * @param groupdField Group of vertices assocated with faces of
   *   cells defining fault surface
   * @param flipFault Flag indicating to flip positive/negative sides
   * of the fault.
   */
  static
  void createFault(topology::SubMesh* faultMesh,
		   ALE::Obj<SieveFlexMesh>& faultBoundary,
		   DM& faultBoundaryDM,
		   const topology::Mesh& mesh,
		   DMLabel groupField,
		   const bool flipFault =false);

  /** Create cohesive cells.
   *
   * If firstFaultVertex == 0, then firstFaultVertex is set to the first point
   * not currently used in the mesh, and firstFaultCell is incremented with this
   * point. These values are updated as new fault vertices and cells are added.
   *
   * @param fault Finite-element mesh of fault (output)
   * @param mesh Finite-element mesh
   * @param materialId Material id for cohesive elements.
   * @param firstFaultVertex The first point eligible to become a new fault vertex
   * @param firstFaultCell The first point eligible to become a new fault cell
   * @param constraintCell True if creating cells constrained with 
   *   Lagrange multipliers that require extra vertices, false otherwise
   */
  static
  void create(topology::Mesh* mesh,
	      const topology::SubMesh& faultMesh,
              const ALE::Obj<SieveFlexMesh>& faultBoundary,
              DM faultBoundaryDM,
              DMLabel groupField,
              const int materialId,
              int& firstFaultVertex,
              int& firstLagrangeVertex,
              int& firstFaultCell,
              const bool constraintCell = false);

  /** Create cohesive cells in an interpolated mesh.
   *
   * If firstFaultVertex == 0, then firstFaultVertex is set to the first point
   * not currently used in the mesh, and firstFaultCell is incremented with this
   * point. These values are updated as new fault vertices and cells are added.
   *
   * @param fault Finite-element mesh of fault (output)
   * @param mesh Finite-element mesh
   * @param materialId Material id for cohesive elements.
   * @param firstFaultVertex The first point eligible to become a new fault vertex
   * @param firstFaultCell The first point eligible to become a new fault cell
   * @param constraintCell True if creating cells constrained with 
   *   Lagrange multipliers that require extra vertices, false otherwise
   */
  static
  void createInterpolated(topology::Mesh* mesh,
                          const topology::SubMesh& faultMesh,
                          const ALE::Obj<SieveFlexMesh>& faultBoundary,
                          const ALE::Obj<topology::Mesh::IntSection>& groupField,
                          const int materialId,
                          int& firstFaultVertex,
                          int& firstLagrangeVertex,
                          int& firstFaultCell,
                          const bool constraintCell = false);

  /** Create (distributed) fault mesh from cohesive cells.
   *
   * @param faultMesh Finite-element mesh of fault (output).
   * @param mesh Finite-element mesh.
   * @param materialId Material id for cohesive elements.
   * @param constraintCell True if creating cells constrained with 
   *   Lagrange multipliers that require extra vertices, false otherwise.
   */
  static
  void createFaultParallel(topology::SubMesh* faultMesh,
			   const topology::Mesh& mesh,
			   const int materialId,
			   const bool constraintCell =false);

}; // class CohesiveTopology

#endif // pylith_faults_cohesivetopology_hh


// End of file 
