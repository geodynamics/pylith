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

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/array.hh" // USES double_array, int_array

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "Selection.hh" // USES boundary()

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::IntSection IntSection;

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIO::MeshIO(void) :
  _mesh(0),
  _debug(false),
  _interpolate(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIO::~MeshIO(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::MeshIO::deallocate(void)
{ // deallocate
} // deallocate
  
// ----------------------------------------------------------------------
// Get spatial dimension of mesh.
int
pylith::meshio::MeshIO::getMeshDim(void) const
{ // getMeshDim
  return (0 != _mesh) ? _mesh->dimension() : 0;
} // getMeshDim

// ----------------------------------------------------------------------
// Read mesh from file.
void 
pylith::meshio::MeshIO::read(topology::Mesh* mesh)
{ // read
  assert(0 == _mesh);

  _mesh = mesh;
  _mesh->debug(_debug);
  _read();

  _mesh = 0;
} // read

// ----------------------------------------------------------------------
// Write mesh to file.
void 
pylith::meshio::MeshIO::write(topology::Mesh* const mesh)
{ // write
  assert(0 == _mesh);

  _mesh = mesh;
  _write();
  _mesh = 0;
} // write

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

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = 
    vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = 
    vertices->end();
  const ALE::Obj<RealSection>& coordsField =
    sieveMesh->hasRealSection("coordinates_dimensioned") ?
    sieveMesh->getRealSection("coordinates_dimensioned") :
    sieveMesh->getRealSection("coordinates");
  assert(!coordsField.isNull());

  *numVertices = vertices->size();
  assert(*numVertices > 0);
  assert(vertices->size() > 0);
  *spaceDim = coordsField->getFiberDimension(*vertices->begin());
  assert(*spaceDim > 0);

  const int size = (*numVertices) * (*spaceDim);
  coordinates->resize(size);

  int i = 0;
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
      v_iter != verticesEnd;
      ++v_iter) {
    const RealSection::value_type *vertexCoords = 
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

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& meshCells = 
    sieveMesh->heightStratum(0);
  assert(!meshCells.isNull());
  const SieveMesh::label_sequence::iterator meshCellsBegin = 
    meshCells->begin();
  const SieveMesh::label_sequence::iterator meshCellsEnd = 
    meshCells->end();

  *meshDim = _mesh->dimension();
  *numCells = meshCells->size();
  *numCorners = sieveMesh->getNumCellCorners();
  
  const ALE::Obj<SieveMesh::numbering_type>& vNumbering = 
    sieveMesh->getFactory()->getLocalNumbering(sieveMesh, 0);

  const int size = (*numCells) * (*numCorners);
  cells->resize(size);
    
  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> 
    pV(sieve->getMaxConeSize());
  int i = 0;
  for(SieveMesh::label_sequence::iterator e_iter = meshCellsBegin;
      e_iter != meshCellsEnd;
      ++e_iter) {
    sieve->cone(*e_iter, pV);
    const SieveMesh::point_type *cone = pV.getPoints();
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

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  assert(!sieveMesh.isNull());

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  ///logger.setDebug(2);
  logger.stagePush("Materials");
  const ALE::Obj<SieveMesh::label_type>& labelMaterials = 
    sieveMesh->createLabel("material-id");
  if (!sieveMesh->commRank()) {
    const ALE::Obj<SieveMesh::label_sequence>& cells = 
      sieveMesh->heightStratum(0);
    assert(!cells.isNull());
    const SieveMesh::label_sequence::iterator cellsBegin = 
      cells->begin();
    const SieveMesh::label_sequence::iterator cellsEnd = 
      cells->end();

    const int numCells = materialIds.size();
    if (cells->size() != numCells) {
      std::ostringstream msg;
      msg << "Mismatch in size of materials identifier array ("
          << numCells << ") and number of cells in mesh ("
          << cells->size() << ").";
      throw std::runtime_error(msg.str());
    } // if
    int i = 0;

#ifdef IMESH_NEW_LABELS
    labelMaterials->setChart(sieveMesh->getSieve()->getChart());
    for(SieveMesh::label_sequence::iterator e_iter = cellsBegin;
	e_iter != cellsEnd;
	++e_iter) {
      labelMaterials->setConeSize(*e_iter, 1);
    }
    if (cells->size()) {labelMaterials->setSupportSize(0, cells->size());}
    labelMaterials->allocate();
#endif
    for(SieveMesh::label_sequence::iterator e_iter = cellsBegin;
	e_iter != cellsEnd;
	++e_iter) {
      sieveMesh->setValue(labelMaterials, *e_iter, materialIds[i++]);
    }
#ifdef IMESH_NEW_LABELS
    labelMaterials->recalculateLabel();
#endif
  } // if
  logger.stagePop();
  ///logger.setDebug(0);
} // _setMaterials

// ----------------------------------------------------------------------
// Get material identifiers for cells.
void
pylith::meshio::MeshIO::_getMaterials(int_array* materialIds) const
{ // _getMaterials
  assert(0 != materialIds);
  assert(0 != _mesh);

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const int numCells = cells->size();

  const int size = numCells;
  materialIds->resize(size);
  
  const ALE::Obj<SieveMesh::label_type>& labelMaterials = 
    sieveMesh->getLabel("material-id");
  const int idDefault = 0;
  
  int i = 0;
  for(SieveMesh::label_sequence::iterator e_iter = cellsBegin;
      e_iter != cellsEnd;
      ++e_iter)
    (*materialIds)[i++] = 
      sieveMesh->getValue(labelMaterials, *e_iter, idDefault);
} // _getMaterials

// ----------------------------------------------------------------------
// Build a point group as an int section.
void
pylith::meshio::MeshIO::_setGroup(const std::string& name,
				  const GroupPtType type,
				  const int_array& points)
{ // _setGroup
  assert(0 != _mesh);

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  assert(!sieveMesh.isNull());

  if (sieveMesh->hasIntSection(name)) {
    std::ostringstream msg;
    msg << "Could not setup group '" << name
	<< "'. Group already exists in mesh.";
    throw std::runtime_error(msg.str());
  } // if

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("VertexGroups");
  const ALE::Obj<IntSection>& groupField = sieveMesh->getIntSection(name);
  assert(!groupField.isNull());

  const int numPoints   = points.size();
  const int numVertices = sieveMesh->depthStratum(0)->size();
  const int numCells    = sieveMesh->heightStratum(0)->size();
  if (CELL == type) {
    groupField->setChart(IntSection::chart_type(0,numCells));
    for(int i=0; i < numPoints; ++i)
      groupField->setFiberDimension(points[i], 1);
  } else if (VERTEX == type) {
    groupField->setChart(IntSection::chart_type(numCells, 
						numCells+numVertices));
    for(int i=0; i < numPoints; ++i)
      groupField->setFiberDimension(numCells+points[i], 1);
  } // if/else
  sieveMesh->allocate(groupField);
  logger.stagePop();
} // _setGroup

// ----------------------------------------------------------------------
// Create empty groups on other processes
void
pylith::meshio::MeshIO::_distributeGroups()
{ // _distributeGroups
  assert(0 != _mesh);

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  assert(!sieveMesh.isNull());

  if (!sieveMesh->commRank()) {
    const ALE::Obj<std::set<std::string> >& sectionNames = 
      sieveMesh->getIntSections();
    int numGroups = sectionNames->size();

    MPI_Bcast(&numGroups, 1, MPI_INT, 0, sieveMesh->comm());

    const std::set<std::string>::const_iterator namesEnd = sectionNames->end();
    for (std::set<std::string>::const_iterator name=sectionNames->begin();
         name != namesEnd;
	 ++name) {
      int len = name->size();
      
      MPI_Bcast(&len, 1, MPI_INT, 0, sieveMesh->comm());
      MPI_Bcast((void *) name->c_str(), len, MPI_CHAR, 0, sieveMesh->comm());
    }
  } else {
    int numGroups;

    MPI_Bcast(&numGroups, 1, MPI_INT, 0, sieveMesh->comm());
    for(int g = 0; g < numGroups; ++g) {
      char *name;
      int   len;

      MPI_Bcast(&len, 1, MPI_INT, 0, sieveMesh->comm());
      name = new char[len+1];
      MPI_Bcast(name, len, MPI_CHAR, 0, sieveMesh->comm());
      name[len] = 0;
      const ALE::Obj<IntSection>& groupField = sieveMesh->getIntSection(name);
      assert(!groupField.isNull());
      sieveMesh->allocate(groupField);
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

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<std::set<std::string> >& sectionNames = 
    sieveMesh->getIntSections();
  
  const int numGroups = sectionNames->size();
  names->resize(numGroups);
  const std::set<std::string>::const_iterator namesEnd = sectionNames->end();
  int i=0;
  for (std::set<std::string>::const_iterator name=sectionNames->begin();
       name != namesEnd;
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

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  assert(!sieveMesh.isNull());

  if (!sieveMesh->hasIntSection(name)) {
    std::ostringstream msg;
    msg << "Could not get group '" << name
	<< "'. Group missing from mesh.";
    throw std::runtime_error(msg.str());
  } // if

  const ALE::Obj<IntSection>& groupField = sieveMesh->getIntSection(name);
  assert(!groupField.isNull());
  const IntSection::chart_type& chart = groupField->getChart();
  SieveMesh::point_type firstPoint;
  IntSection::chart_type::const_iterator chartEnd = chart.end();
  for(IntSection::chart_type::const_iterator c_iter=chart.begin();
      c_iter != chartEnd;
      ++c_iter) {
    if (groupField->getFiberDimension(*c_iter)) {
      firstPoint = *c_iter;
      break;
    } // if
  }
  ALE::Obj<SieveMesh::numbering_type> numbering;

  if (sieveMesh->height(firstPoint) == 0) {
    *type = CELL;
    numbering = sieveMesh->getFactory()->getNumbering(sieveMesh, 
						      sieveMesh->depth());
  } else {
    *type = VERTEX;
    numbering = sieveMesh->getFactory()->getNumbering(sieveMesh, 0);
  } // if/else
  const int numPoints = groupField->size();
  points->resize(numPoints);
  int i = 0;

  for(IntSection::chart_type::const_iterator c_iter=chart.begin();
      c_iter != chartEnd;
      ++c_iter) {
    assert(!numbering.isNull());
    if (groupField->getFiberDimension(*c_iter)) (*points)[i++] = 
      numbering->getIndex(*c_iter);
  } // for
} // _getGroup


// End of file
