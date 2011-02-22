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

#include <portinfo>

#include "UCDFaultFile.hh" // implementation of class methods

#include "MeshBuilder.hh" // USES MeshBuilder

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/utils/array.hh" // USES double_array, int_array

#include <petsc.h> // USES MPI_Comm

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
int
pylith::meshio::UCDFaultFile::numVertices(const char* filename)
{ // read
  int nvertices = 0;

  std::ifstream fin(filename, std::ios::in);
  if (!(fin.is_open() && fin.good())) {
    std::ostringstream msg;
    msg << "Could not open ASCII INP file '" << filename << "' for reading.";
    throw std::runtime_error(msg.str());
  } // if
    
  // Section 1: <num_nodes> <num_cells> <num_ndata> <num_cdata> <num_mdata>
  fin >> nvertices;
  fin.close();

  if (nvertices <= 0) {
    std::ostringstream msg;
    msg << "Number of vertices (" << nvertices << ") in ASCII INP file '"
	<< filename << "' must be positive.";
    throw std::runtime_error(msg.str());
  } // if
  
  return nvertices;
} // numVertices


// ----------------------------------------------------------------------
void
pylith::meshio::UCDFaultFile::read(const char* filename,
				   topology::SubMesh* faultMesh,
				   ALE::Obj<FlexMesh>& faultBoundary,
				   const topology::Mesh& mesh)
{ // read
  assert(0 != faultMesh);

  int faultDim = 2;
  int fSpaceDim = 0;
  int numFVertices = 0;
  int numFCells = 0;
  int numFCorners = 0;
  double_array fCoordinates;
  int_array fCells;
  int_array fMaterialIds;
  int_array faceCells;
  int_array vertexIDs;

  const ALE::Obj<topology::Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  int rank = 0;
  MPI_Comm_rank(mesh.comm(), &rank);
  if (0 == rank) {
    std::ifstream fin(filename, std::ios::in);
    if (!(fin.is_open() && fin.good())) {
      std::ostringstream msg;
      msg << "Could not open ASCII INP file '" << filename << "' for reading.";
      throw std::runtime_error(msg.str());
    } // if
    int numNodeData = 0;
    int numCellData = 0;
    int numModelData = 0;
    
    fSpaceDim = 3;

    // Section 1: <num_nodes> <num_cells> <num_ndata> <num_cdata> <num_mdata>
    fin >> numFVertices;
    fin >> numFCells;
    fin >> numNodeData;
    fin >> numCellData;
    fin >> numModelData;
    // Section 2: <node_id 1> <x> <y> <z>
    fCoordinates.resize(numFVertices * fSpaceDim);
    for(int v = 0; v < numFVertices; ++v) {
      int id;
      fin >> id;
      for(int d = 0; d < fSpaceDim; ++d)
        fin >> fCoordinates[v*fSpaceDim + d];
    } // for

    // Section 3: <cell_id 1> <mat_id> <cell_type> <cell_vert 1> ... <cell_vert n> 
    fMaterialIds.resize(numFCells);
    for(int c = 0; c < numFCells; ++c) {
      std::string cellType = "";
      int cellID = 0;
      
      fin >> cellID;
      fin >> fMaterialIds[c];
      fin >> cellType;
      if (cellType == "tri") {
        numFCorners = 3;
      } else if (cellType == "quad") {
        numFCorners = 4;
      } else {
        std::ostringstream msg;
        msg << "Unknown cell type " << cellType << "while reading INP file.";
        throw std::runtime_error(msg.str());
      } // if/else
      if (c == 0) fCells.resize(numFCells*numFCorners);
      for(int i = 0; i < numFCorners; ++i)
	fin >> fCells[c*numFCorners+i];
    } // for

    // Section 4: <num_comp for node data> <size comp 1> <size comp 2>...<size comp n>
    int numComponents = 0;
    int totalSize = 0;

    fin >> numComponents;
    totalSize = 0;
    for(int i = 0; i < numComponents; ++i) {
      int compSize;
      fin >> compSize;
      totalSize += compSize;
    } // for
    // Section 5: <node_comp_label 1> , <units_label 1>
    for(int c = 0; c < numComponents; ++c) {
      std::string label, typeName;

      fin >> label;
      fin >> typeName;
    } // for

    // Section 6: <node_id 1> <node_data 1> ... <node_data num_ndata>
    vertexIDs.resize(numFVertices);
    for(int v = 0; v < numFVertices; ++v) {
      int id;
      double dummy;

      fin >> id;
      fin >> vertexIDs[v]; // global node number
      fin >> dummy; // components of normal at vertices
      fin >> dummy;
      fin >> dummy;
    } // for

    // Section 7: <num_comp for cell's data> <size comp 1> <size comp 2>...<size comp n>
    fin >> numComponents;
    totalSize = 0;
    for(int i = 0; i < numComponents; ++i) {
      int compSize;
      fin >> compSize;
      totalSize += compSize;
    }

    // Section 8: <cell_component_label 1> , <units_label 1> 
    for(int c = 0; c < numComponents; ++c) {
      std::string label, typeName;

      fin >> label;
      fin >> typeName;
    } // for

    // Section 9: <cell_id 1> <cell_data 1> ... <cell_data num_cdata> 
    faceCells.resize(numFCells*2);
    for(int c = 0; c < numFCells; ++c) {
      int id, faultId;

      fin >> id;
      fin >> faceCells[c*2+0]; // Cell numbers in global mesh on either side of fault
      fin >> faceCells[c*2+1];
      fin >> faultId;
    } // for

    // Determine the number of cells
    //   Only do this once since we add cohesive cells after that
    static int numCells = -1;

    assert(!sieveMesh->heightStratum(0).isNull());
    if (numCells == -1) {numCells = sieveMesh->heightStratum(0)->size();}

    // Renumber vertices and use zero based indices
    // UCD file has one-based indices for both vertexIDs and fCells
    //   Also, vertex numbers are offset by the number of cells
    for(int c = 0; c < numFCells; ++c)
      for(int corner = 0; corner < numFCorners; ++corner)
        fCells[c*numFCorners+corner] = 
	  vertexIDs[fCells[c*numFCorners+corner]-1] - 1 + numCells;
    
    // Switch to zero based index for global cell numbering
    for (int c=0; c < numFCells; ++c)
      for (int i=0; i < 2; ++i)
	faceCells[c*2+i] -= 1;
  } // if

  // Create Sieve mesh for fault
  ALE::Obj<pylith::topology::Mesh::SieveSubMesh>& faultSieveMesh =
    faultMesh->sieveMesh();
  faultSieveMesh =
    new pylith::topology::Mesh::SieveSubMesh(mesh.comm(), mesh.dimension()-1,
					     mesh.debug());
  
  assert(!sieveMesh->getSieve().isNull());
  const int firstFaultCell = 
    sieveMesh->getSieve()->getBaseSize() + sieveMesh->getSieve()->getCapSize();
  MeshBuilder::buildFaultMesh(faultSieveMesh, faultBoundary, 
			      fCoordinates, numFVertices, fSpaceDim, fCells, 
			      numFCells, numFCorners, firstFaultCell, 
			      faceCells, faultDim);
} // read


// End of file

