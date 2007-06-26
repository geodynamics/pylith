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

#include "MeshIO.hh" // MeshIOHDF5 ISA MeshIO
#include "MeshIOHDF5.hh" // implementation of class methods

#include "PetscMesh.hh"
extern "C" {
#include "hdf5.h" 
}

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
MeshIOHDF5::MeshIOHDF5(void) :
  _filename("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
MeshIOHDF5::~MeshIOHDF5(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Read mesh from file.
void
MeshIOHDF5::read(ALE::Obj<ALE::PetscMesh>* pMesh)
{ // read
  assert(0 != pMesh);
  MPI_Comm comm = PETSC_COMM_WORLD;
  int rank;
  int meshDim = 0;
  int numDims = 0;
  int numVertices = 0;
  int numElements = 0;
  int numCorners = 0;
  double* coordinates = 0;
  int* elements = 0;

  MPI_Comm_rank(comm, &rank);
  if (!rank) {
    hid_t filein = H5Fopen(_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (filein < 0) {
      std::ostringstream msg;
      msg << "Could not open HDF5 mesh file '" << _filename
          << "' for reading.\n";
      throw std::runtime_error(msg.str());
    } // if

    _readMeshInfo(filein, pMesh);
    _readVertices(filein, &coordinates, &numVertices, &numDims);
    _readElements(filein, &elements, &numElements, &numCorners);
  
    *pMesh = ALE::PetscMesh(PETSC_COMM_WORLD, meshDim);
    (*pMesh)->debug = true;
    bool interpolate = false;

    _buildMesh(coordinates, numVertices, spaceDim,
               cells, numCells, numCorners, meshDim);
    delete[] coordinates; coordinates = 0;
    delete[] elements; elements = 0;

    // loop over charts
    // _readChart();

    H5Fclose(filein);
  } else {
    _buildMesh(coordinates, numVertices, spaceDim,
               cells, numCells, numCorners, meshDim);
  }
  //_setMaterials(materialIds);
  //_distributeGroups();
} // read

// ----------------------------------------------------------------------
// Write mesh to file.
void
MeshIOHDF5::write(const ALE::Obj<ALE::PetscMesh>& mesh) const
{ // write
  HDF5 fileout(_filename.c_str(), H5F_ACC_TRUNC);

  try {
    _writeMeshInfo(fileout, mesh);
    _writeVertices(fileout, mesh);
    _writeElements(fileout, mesh);

    // LOOP OVER CHARTS
    // _writeChart(fileout, mesh, nameIter->c_str());

  } catch (std::exception& err) {
    std::ostringstream msg;
    msg
      << "Error occurred while writing HDF5 file '" << _filename << "'.\n"
      << err.what();
    throw std::runtime_error(msg.str());
  } catch (...) {
    std::ostringstream msg;
    msg
      << "Unknown error occurred while writing HDF5 file '" 
      << _filename << "'";
    throw std::runtime_error(msg.str());
  } // try/catch
} // write

// ----------------------------------------------------------------------
// Read general mesh information.
void
MeshIOHDF5::_readMeshInfo(hid_t& filein,
			  ALE::Obj<ALE::PetscMesh>* pMesh)
{ // _readMeshInfo
} // _readMeshInfo

// ----------------------------------------------------------------------
// Write general mesh information.
void
MeshIOHDF5::_writeMeshInfo(HDF5& fileout,
			   const ALE::Obj<ALE::PetscMesh>& mesh) const
{ // _writeMeshInfo
  hid_t meshGroup = fileout.createGroup("/mesh");

  const int dimension = mesh->getDimension();
  fileout.createAttribute(meshGroup, "dimension", 
			  (void*) &dimension, H5T_NATIVE_INT);

  if (H5Gclose(meshGroup) < 0)
    throw std::runtime_error("Could not close 'mesh' group. ");
} // _writeMeshInfo

// ----------------------------------------------------------------------
// Read mesh vertices.
void
MeshIOHDF5::_readVertices(hid_t& filein,
			  double** pCoordinates,
			  int* pNumVertices, 
			  int* pNumDims) const
{ // _readVertices
  double* coordinates = 0;
  int numDims = 0;
  int numVertices = 0;

  // ADD STUFF HERE

  if (0 != pCoordinates)
    *pCoordinates = coordinates;
  if (0 != pNumVertices)
    *pNumVertices = numVertices;
  if (0 != pNumDims)
    *pNumDims = numDims;
} // _readVertices

// ----------------------------------------------------------------------
// Write mesh vertices.
void
MeshIOHDF5::_writeVertices(HDF5& fileout,
			   const ALE::Obj<ALE::PetscMesh>& mesh) const
{ // _writeVertices
  ALE::Obj<ALE::PetscMesh::field_type> coords_field = mesh->getCoordinates();
  ALE::Obj<ALE::PetscMesh::bundle_type> vertexBundle = mesh->getBundle(0);
  ALE::PetscMesh::field_type::patch_type patch;
  const double* coordinates = coords_field->restrict(patch);
  const int numVertices = (vertexBundle->getGlobalOffsets()) ?
    vertexBundle->getGlobalOffsets()[mesh->commSize()] :
    mesh->getTopology()->depthStratum(0)->size();
  const int numDims = coords_field->getFiberDimension(patch, 
			      *mesh->getTopology()->depthStratum(0)->begin());

  hid_t verticesId = fileout.createGroup("/mesh/vertices");
  int dims[2];
  dims[0] = numVertices;
  dims[1] = numDims;
  fileout.writeDataset(verticesId, "coordinates", (void*) coordinates, 
		       dims, 2, H5T_NATIVE_DOUBLE);
  
  if (H5Gclose(verticesId) < 0)
    throw std::runtime_error("Could not close 'mesh/vertices' group. ");
} // _writeVertices
  
// ----------------------------------------------------------------------
// Read mesh elements.
void
MeshIOHDF5::_readElements(hid_t& filein,
			   int** pElements,
			   int* pNumElements, 
			   int* pNumCorners) const
{ // _readElements
  int* elements = 0;
  int numElements = 0;
  int numCorners = 0;
  int dimension = 0;

  // ADD STUFF HERE

  if (0 != pElements)
    *pElements = elements;
  if (0 != pNumElements)
    *pNumElements = numElements;
  if (0 != pNumCorners)
    *pNumCorners = numCorners;
} // _readElements

// ----------------------------------------------------------------------
// Write mesh elements.
void
MeshIOHDF5::_writeElements(hid_t& fileout,
			   const ALE::Obj<ALE::PetscMesh>& mesh) const
{ // _writeElements
  ALE::Obj<ALE::PetscMesh::sieve_type> topology = mesh->getTopology();
  ALE::Obj<ALE::PetscMesh::sieve_type::traits::heightSequence> elements = 
    topology->heightStratum(0);
  ALE::Obj<ALE::PetscMesh::bundle_type> vertexBundle =  mesh->getBundle(0);
  ALE::PetscMesh::bundle_type::patch_type patch;
  std::string orderName("element");

  assert(0 != topology);
  assert(0 != elements);
  assert(0 != vertexBundle);

  int numCorners = 
    topology->nCone(*elements->begin(), topology->depth())->size();
  const int numElements = mesh->getTopology()->heightStratum(0)->size();

  int* simplices = 
    (numElements*numCorners > 0) ? new int[numElements*numCorners] : 0;
  const int offset = (useIndexZero()) ? 0 : 1;
  int index = 0;
  for(ALE::PetscMesh::sieve_type::traits::heightSequence::iterator e_itor = 
	elements->begin(); 
      e_itor != elements->end();
      ++e_itor) {
    ALE::Obj<ALE::PetscMesh::bundle_type::order_type::coneSequence> cone = 
      vertexBundle->getPatch(orderName, *e_itor);
    assert(0 != cone);
    for(ALE::PetscMesh::bundle_type::order_type::coneSequence::iterator c_itor = 
	  cone->begin(); 
	c_itor != cone->end(); 
	++c_itor)
      simplices[index++] = 
	offset + vertexBundle->getIndex(patch, *c_itor).prefix;

  hid_t elementsId = fileout.createGroup("/mesh/elements");
  int dims[2];
  dims[0] = numElements;
  dims[1] = numCorners;
  fileout.writeDataset(elementsId, "simplices", (void*) simplices, 
		       dims, 2, H5T_NATIVE_INT);
  
  if (H5Gclose(elementsId) < 0)
    throw std::runtime_error("Could not close 'mesh/vertices' group. ");
} // _writeElements

// ----------------------------------------------------------------------
// Read mesh charts.
void
MeshIOHDF5::_readChart(hid_t& filein,
		       ALE::Obj<ALE::PetscMesh>* pMesh) const
{ // _readChart
  std::string name = ""; // Name of chart
  int dimension = 0; // Topology dimension associated with chart
  int count = 0; // Number of entities in chart
  int* indices = 0; // Indices of entities in chart

} // _readChart

// ----------------------------------------------------------------------
// Write mesh chart.
void
MeshIOHDF5::_writeChart(hid_t& fileout,
			 const ALE::Obj<ALE::PetscMesh>& mesh,
			 const char* name) const
{ // _writeChart
} // _writeChart
  
// End of file 
