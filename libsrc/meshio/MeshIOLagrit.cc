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

#include "MeshIOLagrit.hh" // implementation of class methods

#include "GMVFileAscii.hh" // USES GMVFileAscii
#include "GMVFileBinary.hh" // USES GMVFileBinary
#include "PsetFileAscii.hh" // USES PsetFileAscii
#include "PsetFileBinary.hh" // USES PsetFileBinary

#include "pylith/utils/array.hh" // USES double_array, int_array

#include <petsc.h> // USES MPI_Comm

#include <assert.h> // USES assert()
#include <stdexcept> // TEMPORARY

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOLagrit::MeshIOLagrit(void) :
  _filenameGmv(""),
  _filenamePset(""),
  _flipEndian(false),
  _ioInt32(true),
  _isRecordHeader32Bit(true)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIOLagrit::~MeshIOLagrit(void)
{ // destructor
} // destructor

void
pylith::meshio::MeshIOLagrit::readFault(const std::string filename, const Obj<Mesh>& mesh, const Obj<Mesh>& fault, Obj<ALE::Mesh>& faultBd) {
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
    std::ifstream fin(filename.c_str(), std::ios::in);
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
      for(int d = 0; d < fSpaceDim; ++d) {
        fin >> fCoordinates[v*fSpaceDim + d];
      }
    }
    // Section 3: <cell_id 1> <mat_id> <cell_type> <cell_vert 1> ... <cell_vert n> 
    fMaterialIds.resize(numFCells);
    for(int c = 0; c < numFCells; ++c) {
      std::string cellType;
      int         cellID;

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
      }
      if (c == 0) fCells.resize(numFCells*numFCorners);
      for(int i = 0; i < numFCorners; ++i) fin >> fCells[c*numFCorners+i];
    }
    // Section 4: <num_comp for node data> <size comp 1> <size comp 2>...<size comp n>
    int numComponents, totalSize;

    fin >> numComponents;
    totalSize = 0;
    for(int i = 0; i < numComponents; ++i) {
      int compSize;
      fin >> compSize;
      totalSize += compSize;
    }
    // Section 5: <node_comp_label 1> , <units_label 1>
    for(int c = 0; c < numComponents; ++c) {
      std::string label, typeName;

      fin >> label;
      fin >> typeName;
    }
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
    }
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
    }
    // Section 9: <cell_id 1> <cell_data 1> ... <cell_data num_cdata> 
    faceCells.resize(numFCells*2);
    for(int c = 0; c < numFCells; ++c) {
      int id, dummy;

      fin >> id;
      fin >> faceCells[c*2+0]; // Cell numbers in global mesh on either side of fault
      fin >> faceCells[c*2+1];
    }

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
  }

  const int firstFaultCell = mesh->getSieve()->getBaseSize() + mesh->getSieve()->getCapSize();
  _buildFaultMesh(fCoordinates, numFVertices, fSpaceDim, fCells, numFCells, numFCorners, firstFaultCell, faceCells, faultDim, fault, faultBd);
}

// ----------------------------------------------------------------------
// Unpickle mesh
void
pylith::meshio::MeshIOLagrit::_read(void)
{ // _read
  MPI_Comm comm = PETSC_COMM_WORLD;
  int rank;
  int meshDim = 0;
  int spaceDim = 0;
  int numVertices = 0;
  int numCells = 0;
  int numCorners = 0;
  double_array coordinates;
  int_array cells;
  int_array materialIds;

  MPI_Comm_rank(comm, &rank);
  if (!rank) {
    if (GMVFile::isAscii(_filenameGmv.c_str())) {
      GMVFileAscii filein(_filenameGmv.c_str());
      filein.read(&coordinates, &cells, &materialIds, 
                  &meshDim, &spaceDim, &numVertices, &numCells, &numCorners);
      _orientCellsAscii(&cells, numCells, numCorners, meshDim);
    } else {
      GMVFileBinary filein(_filenameGmv.c_str(), _flipEndian);
      filein.read(&coordinates, &cells, &materialIds, 
                  &meshDim, &spaceDim, &numVertices, &numCells, &numCorners);
      _orientCellsBinary(&cells, numCells, numCorners, meshDim);
    } // if/else
  }
  _buildMesh(coordinates, numVertices, spaceDim,
             cells, numCells, numCorners, meshDim);
  _setMaterials(materialIds);

  if (!rank) {
    std::vector<PsetFile::Pset> groups;
    if (PsetFile::isAscii(_filenamePset.c_str())) {
      PsetFileAscii filein(_filenamePset.c_str());
      filein.read(&groups);
    } else {
      PsetFileBinary filein(_filenamePset.c_str(), 
			    _flipEndian,
			    _ioInt32,
			    _isRecordHeader32Bit);
      filein.read(&groups);
    } // if/else
    GroupPtType type = VERTEX;
    const int numGroups = groups.size();
    for (int iGroup=0; iGroup < numGroups; ++iGroup)
      _setGroup(groups[iGroup].name, type, groups[iGroup].points);
  }
  _distributeGroups();
} // _read

// ----------------------------------------------------------------------
// Pickle mesh
void
pylith::meshio::MeshIOLagrit::_write(void) const
{ // _write
  throw std::logic_error("MeshIOLagrit::_write not implemented.");
} // _write

// ----------------------------------------------------------------------
// Reorder vertices in cells from ASCII GMV file to match PyLith
// conventions.
void
pylith::meshio::MeshIOLagrit::_orientCellsAscii(int_array* const cells,
						const int numCells,
						const int numCorners,
						const int meshDim)
{ // _orientCellsAscii
  assert(0 != cells);
  assert(cells->size() == numCells*numCorners);

  if (3 == meshDim && 4 == numCorners) // TET
    for (int iCell=0; iCell < numCells; ++iCell) {
      const int i1 = iCell*numCorners+1;
      const int i2 = iCell*numCorners+2;
      const int tmp = (*cells)[i1];
      (*cells)[i1] = (*cells)[i2];
      (*cells)[i2] = tmp;
    } // for
} // _orientCellsAscii
  
// ----------------------------------------------------------------------
// Reorder vertices in cells from binary GMV file to match PyLith
// conventions.
void
pylith::meshio::MeshIOLagrit::_orientCellsBinary(int_array* const cells,
						 const int numCells,
						 const int numCorners,
						 const int meshDim)
{ // _orientCellsBinary
  assert(0 != cells);
  assert(cells->size() == numCells*numCorners);

  if (3 == meshDim && 4 == numCorners)  // TET
    ; // do nothing
} // _orientCellsBinary
  

// End of file 
