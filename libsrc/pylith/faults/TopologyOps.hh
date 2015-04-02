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
// Copyright (c) 2010-2015 University of California, Davis
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
#include <set>
#include <vector>
#include <iostream>

// TopologyOps ----------------------------------------------------------
/// Helper object for creation of cohesive cells.
class pylith::faults::TopologyOps
{ // class TopologyOps

  // PUBLIC TYPEDEFS ////////////////////////////////////////////////////
public :

  typedef std::set<PetscInt> PointSet;
  typedef std::vector<PetscInt> PointArray;

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  static
  void classifyCellsDM(PetscDM dmMesh,
		       PetscInt vertex,
		       const int depth,
		       const int faceSize,
		       PetscInt firstCohesiveCell,
		       PointSet& replaceCells,
		       PointSet& noReplaceCells,
		       const int debug);

}; // class CohesiveTopology

#endif // pylith_faults_cohesivetopology_hh


// End of file 
