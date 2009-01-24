// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "UCDFaultFile.hh" // implementation of class methods

#include "MeshBuilder.hh" // USES MeshBuilder

#include "pylith/utils/array.hh" // USES double_array, int_array

#include <petsc.h> // USES MPI_Comm

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
void
pylith::meshio::UCDFaultFile::read(const char* filename,
				   const ALE::Obj<SieveMesh>& mesh, 
				   const ALE::Obj<SieveMesh>& fault, 
				   ALE::Obj<ALE::Mesh>& faultBd)
{ // read
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
  
  if (!fault->commRank()) {
    std::ifstream fin(filename, std::ios::in);
    if (!(fin.is_open() && fin.good())) {
      std::ostringstream msg;
      msg << "Could not open ASCII INP file '" << filename << "' for reading.";
      throw std::runtime_error(msg.str());
    } // if
    int numNodeData, numCellData, numModelData;

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
    int numComponents, totalSize;

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

    if (numCells == -1) {numCells = mesh->heightStratum(0)->size();}

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

  const int firstFaultCell = 
    mesh->getSieve()->getBaseSize() + mesh->getSieve()->getCapSize();
  MeshBuilder::buildFaultMesh(fault, faultBd, 
			      fCoordinates, numFVertices, fSpaceDim, fCells, 
			      numFCells, numFCorners, firstFaultCell, 
			      faceCells, faultDim);
} // read


// End of file

