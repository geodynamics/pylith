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

#include "MeshIO.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES double_array, int_array

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIO::MeshIO(void) :
  _debug(false),
  _interpolate(false),
  _mesh(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIO::~MeshIO(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Get spatial dimension of mesh.
int
pylith::meshio::MeshIO::getMeshDim(void) const
{ // getMeshDim
  assert(0 != _mesh);
  assert(!_mesh->isNull());
  return (*_mesh)->getDimension();
} // getMeshDim

// ----------------------------------------------------------------------
// Read mesh from file.
void 
pylith::meshio::MeshIO::read(ALE::Obj<Mesh>* mesh)
{ // read
  assert(0 == _mesh);

  _mesh = mesh;
  _read();
  _mesh = 0;
  (*mesh)->getFactory()->clear();
} // read

// ----------------------------------------------------------------------
// Write mesh to file.
void 
pylith::meshio::MeshIO::write(ALE::Obj<Mesh>* mesh)
{ // write
  assert(0 == _mesh);

  _mesh = mesh;
  _write();
  _mesh = 0;
} // write

// ----------------------------------------------------------------------
// Set vertices in mesh.
void
pylith::meshio::MeshIO::_buildMesh(const double_array& coordinates,
				   const int numVertices,
				   const int spaceDim,
				   const int_array& cells,
				   const int numCells,
				   const int numCorners,
				   const int meshDim)
{ // _buildMesh
  assert(0 != _mesh);
  MPI_Comm comm = PETSC_COMM_WORLD;
  int      dim  = meshDim;
  int      rank;

  MPI_Bcast(&dim, 1, MPI_INT, 0, comm);
  // :BUG: This causes a memory leak.
  *_mesh = new Mesh(PETSC_COMM_WORLD, dim);

  assert(!_mesh->isNull());
  (*_mesh)->setDebug(_debug);

  ALE::Obj<sieve_type> sieve = new sieve_type((*_mesh)->comm());
  (*_mesh)->setSieve(sieve);

  MPI_Comm_rank(comm, &rank);
  // Memory debugging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.setDebug(_debug);

  logger.stagePush("MeshCreation");
  if (!rank) {
    assert(coordinates.size() == numVertices*spaceDim);
    assert(cells.size() == numCells*numCorners);
    //if (!_interpolate) {
    if (0) {
      // Create the ISieve
      sieve->setChart(Mesh::sieve_type::chart_type(0, numCells+numVertices));
      // Set cone and support sizes
      for(int c = 0; c < numCells; ++c) {sieve->setConeSize(c, numCorners);}
      sieve->symmetrizeSizes(numCells, numCorners, const_cast<int*>(&cells[0]));
      // Allocate point storage
      sieve->allocate();
      // Fill up cones
      int *cone = new int[numCorners];
      for(int c = 0; c < numCells; ++c) {
        for(int v = 0; v < numCorners; ++v) cone[v] = cells[c*numCorners+v]+numCells;
        sieve->setCone(cone, c);
      }
      delete [] cone;
      // Symmetrize to fill up supports
      sieve->symmetrize();
    } else {
      // Same old thing
      ALE::Obj<ALE::Mesh::sieve_type> s = new ALE::Mesh::sieve_type(sieve->comm(), sieve->debug());

      ALE::SieveBuilder<ALE::Mesh>::buildTopology(s, meshDim, 
                                                  numCells, 
                                                  const_cast<int*>(&cells[0]), 
                                                  numVertices, 
                                                  _interpolate,
                                                  numCorners);
      std::map<Mesh::point_type,Mesh::point_type> renumbering;
      ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
    }
    logger.stagePop();
    logger.stagePush("MeshStratification");
    if (!_interpolate) {
      // Optimized stratification
      const ALE::Obj<Mesh::label_type>& height = (*_mesh)->createLabel("height");
      const ALE::Obj<Mesh::label_type>& depth  = (*_mesh)->createLabel("depth");
      for(int c = 0; c < numCells; ++c) {
        height->setCone(0, c);
        depth->setCone(1, c);
      }
      for(int v = numCells; v < numCells+numVertices; ++v) {
        height->setCone(1, v);
        depth->setCone(0, v);
      }
      (*_mesh)->setHeight(1);
      (*_mesh)->setDepth(1);
    } else {
      (*_mesh)->stratify();
    }
    logger.stagePop();
  } else {
    logger.stagePush("MeshStratification");
    (*_mesh)->getSieve()->setChart(sieve_type::chart_type());
    (*_mesh)->getSieve()->allocate();
    (*_mesh)->stratify();
    logger.stagePop();
  }

#if defined(ALE_MEM_LOGGING)
  std::cout
    << std::endl
    << "MeshCreation " << logger.getNumAllocations("MeshCreation")
    << " allocations " << logger.getAllocationTotal("MeshCreation")
    << " bytes" << std::endl
    
    << "MeshCreation " << logger.getNumDeallocations("MeshCreation")
    << " deallocations " << logger.getDeallocationTotal("MeshCreation")
    << " bytes" << std::endl
    
    << "MeshStratification " << logger.getNumAllocations("MeshStratification")
    << " allocations " << logger.getAllocationTotal("MeshStratification")
    << " bytes" << std::endl
    
    << "MeshStratification " << logger.getNumDeallocations("MeshStratification")
    << " deallocations " << logger.getDeallocationTotal("MeshStratification")
    << " bytes" << std::endl << std::endl;
#endif

  ALE::SieveBuilder<Mesh>::buildCoordinates(*_mesh, spaceDim, &coordinates[0]);
} // _buildMesh

// ----------------------------------------------------------------------
// Get coordinates of vertices in mesh.
void
pylith::meshio::MeshIO::_getVertices(double_array* coordinates,
				     int* numVertices,
				     int* spaceDim) const
{ // _getVertices
  assert(0 != coordinates);
  assert(0 != numVertices);
  assert(0 != spaceDim);
  assert(0 != _mesh);
  assert(!_mesh->isNull());

  const ALE::Obj<Mesh::label_sequence>& vertices = (*_mesh)->depthStratum(0);
  assert(!vertices.isNull());
  const ALE::Obj<real_section_type>& coordsField =
    (*_mesh)->getRealSection("coordinates");
  assert(!coordsField.isNull());

  *numVertices = vertices->size();
  assert(*numVertices > 0);
  *spaceDim = coordsField->getFiberDimension(*vertices->begin());
  assert(*spaceDim > 0);

  const int size = (*numVertices) * (*spaceDim);
  coordinates->resize(size);

  int i = 0;
  for(Mesh::label_sequence::iterator v_iter = vertices->begin();
      v_iter != vertices->end();
      ++v_iter) {
    const real_section_type::value_type *vertexCoords = 
      coordsField->restrictPoint(*v_iter);
    for (int iDim=0; iDim < *spaceDim; ++iDim)
      (*coordinates)[i++] = vertexCoords[iDim];
  } // for
} // _getVertices

// ----------------------------------------------------------------------
// Get cells in mesh.
void
pylith::meshio::MeshIO::_getCells(int_array* cells,
				  int* numCells,
				  int* numCorners,
				  int* meshDim) const
{ // _getCells
  assert(0 != cells);
  assert(0 != numCells);
  assert(0 != meshDim);
  assert(0 != _mesh);
  assert(!_mesh->isNull());

  const ALE::Obj<sieve_type>& sieve = (*_mesh)->getSieve();
  assert(!sieve.isNull());
  const ALE::Obj<Mesh::label_sequence>& meshCells = (*_mesh)->heightStratum(0);
  assert(!meshCells.isNull());

  *meshDim = (*_mesh)->getDimension();
  *numCells = meshCells->size();
  *numCorners = (*_mesh)->getNumCellCorners();
  
  const ALE::Obj<Mesh::numbering_type>& vNumbering = 
    (*_mesh)->getFactory()->getLocalNumbering(*_mesh, 0);

  const int size = (*numCells) * (*numCorners);
  cells->resize(size);
    
  ALE::ISieveVisitor::PointRetriever<Mesh::sieve_type> pV(sieve->getMaxConeSize());
  int i = 0;
  for(Mesh::label_sequence::iterator e_iter = meshCells->begin();
      e_iter != meshCells->end();
      ++e_iter) {
    sieve->cone(*e_iter, pV);
    const Mesh::point_type *cone = pV.getPoints();
    for(int p = 0; p < pV.getSize(); ++p, ++i) {
      (*cells)[i] = vNumbering->getIndex(cone[p]);
    }
    pV.clear();
  } // for
} // _getCells

// ----------------------------------------------------------------------
// Tag cells in mesh with material identifiers.
void
pylith::meshio::MeshIO::_setMaterials(const int_array& materialIds)
{ // _setMaterials
  assert(0 != _mesh);
  assert(!_mesh->isNull());

  const ALE::Obj<Mesh::label_type>& labelMaterials = 
    (*_mesh)->createLabel("material-id");
  if (!(*_mesh)->commRank()) {
    const ALE::Obj<Mesh::label_sequence>& cells = (*_mesh)->heightStratum(0);
    assert(!cells.isNull());

    const int numCells = materialIds.size();
    if (cells->size() != numCells) {
      std::ostringstream msg;
      msg << "Mismatch in size of materials identifier array ("
          << numCells << ") and number of cells in mesh ("
          << cells->size() << ").";
      throw std::runtime_error(msg.str());
    } // if
    int i = 0;
    const Mesh::label_sequence::iterator end = cells->end();
    for(Mesh::label_sequence::iterator e_iter = cells->begin();
        e_iter != end;
        ++e_iter)
      (*_mesh)->setValue(labelMaterials, *e_iter, materialIds[i++]);
  } // if
} // _setMaterials

// ----------------------------------------------------------------------
// Get material identifiers for cells.
void
pylith::meshio::MeshIO::_getMaterials(int_array* materialIds) const
{ // _getMaterials
  assert(0 != materialIds);
  assert(0 != _mesh);
  assert(!_mesh->isNull());

  const ALE::Obj<Mesh::label_sequence>& cells = (*_mesh)->heightStratum(0);
  assert(!cells.isNull());
  const int numCells = cells->size();

  const int size = numCells;
  materialIds->resize(size);
  
  const ALE::Obj<Mesh::label_type>& labelMaterials = 
    (*_mesh)->getLabel("material-id");
  const int idDefault = 0;

  int i = 0;
  for(Mesh::label_sequence::iterator e_iter = cells->begin();
      e_iter != cells->end();
      ++e_iter)
    (*materialIds)[i++] = 
      (*_mesh)->getValue(labelMaterials, *e_iter, idDefault);
} // _getMaterials

// ----------------------------------------------------------------------
// Build a point group as an int section.
void
pylith::meshio::MeshIO::_setGroup(const std::string& name,
				  const GroupPtType type,
				  const int_array& points)
{ // _setGroup
  assert(0 != _mesh);
  assert(!_mesh->isNull());

  const ALE::Obj<int_section_type>& groupField = (*_mesh)->getIntSection(name);
  assert(!groupField.isNull());

  const int numPoints   = points.size();
  const int numVertices = (*_mesh)->depthStratum(0)->size();
  const int numCells    = (*_mesh)->heightStratum(0)->size();
  if (CELL == type) {
    groupField->setChart(int_section_type::chart_type(0,numCells));
    for(int i=0; i < numPoints; ++i)
      groupField->setFiberDimension(points[i], 1);
  } else if (VERTEX == type) {
    groupField->setChart(int_section_type::chart_type(numCells,numCells+numVertices));
    for(int i=0; i < numPoints; ++i)
      groupField->setFiberDimension(numCells+points[i], 1);
  } // if/else
  (*_mesh)->allocate(groupField);
} // _setGroup

// ----------------------------------------------------------------------
// Create empty groups on other processes
void
pylith::meshio::MeshIO::_distributeGroups()
{ // _distributeGroups
  assert(0 != _mesh);
  assert(!_mesh->isNull());

  if (!(*_mesh)->commRank()) {
    const ALE::Obj<std::set<std::string> >& sectionNames = 
      (*_mesh)->getIntSections();
    int numGroups = sectionNames->size();

    MPI_Bcast(&numGroups, 1, MPI_INT, 0, (*_mesh)->comm());
    for (std::set<std::string>::const_iterator name=sectionNames->begin();
         name != sectionNames->end(); ++name) {
      int len = name->size();
      
      MPI_Bcast(&len, 1, MPI_INT, 0, (*_mesh)->comm());
      MPI_Bcast((void *) name->c_str(), len, MPI_CHAR, 0, (*_mesh)->comm());
    }
  } else {
    int numGroups;

    MPI_Bcast(&numGroups, 1, MPI_INT, 0, (*_mesh)->comm());
    for(int g = 0; g < numGroups; ++g) {
      char *name;
      int   len;

      MPI_Bcast(&len, 1, MPI_INT, 0, (*_mesh)->comm());
      name = new char[len+1];
      MPI_Bcast(name, len, MPI_CHAR, 0, (*_mesh)->comm());
      name[len] = 0;
      const ALE::Obj<int_section_type>& groupField = (*_mesh)->getIntSection(name);
      assert(!groupField.isNull());
      (*_mesh)->allocate(groupField);
      delete [] name;
    } // for
  } // if/else
} // _distributeGroups

// ----------------------------------------------------------------------
// Get names of all groups in mesh.
void
pylith::meshio::MeshIO::_getGroupNames(string_vector* names) const
{ // _getGroups
  assert(0 != names);
  assert(0 != _mesh);
  assert(!_mesh->isNull());
  
  const ALE::Obj<std::set<std::string> >& sectionNames = 
    (*_mesh)->getIntSections();
  
  const int numGroups = sectionNames->size();
  names->resize(numGroups);
  int i=0;
  for (std::set<std::string>::const_iterator name=sectionNames->begin();
       name != sectionNames->end();
       ++name)
    (*names)[i++] = *name;
} // _getGroups

// ----------------------------------------------------------------------
// Get group entities
void
pylith::meshio::MeshIO::_getGroup(int_array* points,
				  GroupPtType* type,
				  const char *name) const
{ // _getGroup
  assert(0 != points);
  assert(0 != type);
  assert(0 != _mesh);
  assert(!_mesh->isNull());

  const ALE::Obj<int_section_type>& groupField = (*_mesh)->getIntSection(name);
  assert(!groupField.isNull());
  const int_section_type::chart_type& chart = groupField->getChart();
  Mesh::point_type firstPoint;
  for(int_section_type::chart_type::const_iterator c_iter = chart.begin(); c_iter != chart.end(); ++c_iter) {
    if (groupField->getFiberDimension(*c_iter)) {firstPoint = *c_iter; break;}
  }
  ALE::Obj<Mesh::numbering_type> numbering;

  if ((*_mesh)->height(firstPoint) == 0) {
    *type = CELL;
    numbering = (*_mesh)->getFactory()->getNumbering(*_mesh, 
                                                    (*_mesh)->depth());
  } else {
    *type = VERTEX;
    numbering = (*_mesh)->getFactory()->getNumbering(*_mesh, 0);
  } // if/else
  const int numPoints = groupField->size();
  points->resize(numPoints);
  int i = 0;

  for(int_section_type::chart_type::const_iterator c_iter = chart.begin();
      c_iter != chart.end();
      ++c_iter) {
    assert(!numbering.isNull());
    if (groupField->getFiberDimension(*c_iter)) (*points)[i++] = numbering->getIndex(*c_iter);
  } // for
} // _getGroup


// End of file 
