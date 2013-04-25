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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "MeshIO.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "pylith/utils/array.hh" // USES scalar_array, int_array
#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

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
  return (_mesh) ? _mesh->dimension() : 0;
} // getMeshDim

// ----------------------------------------------------------------------
// Read mesh from file.
void 
pylith::meshio::MeshIO::read(topology::Mesh* mesh)
{ // read
  PYLITH_METHOD_BEGIN;

  assert(mesh);
  assert(!_mesh);

  _mesh = mesh;
  _mesh->debug(_debug);
  _read();

  _mesh = 0;

  PYLITH_METHOD_END;
} // read

// ----------------------------------------------------------------------
// Write mesh to file.
void 
pylith::meshio::MeshIO::write(topology::Mesh* const mesh)
{ // write
  PYLITH_METHOD_BEGIN;

  assert(mesh);
  assert(!_mesh);

  _mesh = mesh;
  _write();
  _mesh = 0;

  PYLITH_METHOD_END;
} // write

// ----------------------------------------------------------------------
// Get coordinates of vertices in mesh.
void
pylith::meshio::MeshIO::_getVertices(scalar_array* coordinates,
				     int* numVertices,
				     int* spaceDim) const
{ // _getVertices
  PYLITH_METHOD_BEGIN;

  assert(coordinates);
  assert(numVertices);
  assert(spaceDim);
  assert(_mesh);

#if 0 // Sieve
  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<SieveMesh::label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
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
#endif

  const spatialdata::geocoords::CoordSys* cs = _mesh->coordsys();assert(cs);
  *spaceDim = cs->spaceDim();

  PetscDM dmMesh = _mesh->dmMesh();assert(dmMesh);
  PetscVec coordVec = NULL;
  PetscScalar* coordArray = NULL;
  PetscInt coordSize = 0;    
  PylithScalar lengthScale = 1.0;
  PetscErrorCode err = 0;

  // Get length scale for dimensioning
  err = DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale);

  // Get coordinates and dimensionalize values
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);PYLITH_CHECK_ERROR(err);
  err = VecGetArray(coordVec, &coordArray);PYLITH_CHECK_ERROR(err);
  err = VecGetLocalSize(coordVec, &coordSize);PYLITH_CHECK_ERROR(err);
  assert(coordSize % *spaceDim == 0);
  *numVertices = coordSize / *spaceDim;

  coordinates->resize(coordSize);
  for (PetscInt i=0; i < coordSize; ++i) {
    (*coordinates)[i] = coordArray[i]*lengthScale;
  } // for
  err = VecRestoreArray(coordVec, &coordArray);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _getVertices

// ----------------------------------------------------------------------
// Get cells in mesh.
void
pylith::meshio::MeshIO::_getCells(int_array* cells,
				  int* numCells,
				  int* numCorners,
				  int* meshDim) const
{ // _getCells
  PYLITH_METHOD_BEGIN;

  assert(cells);
  assert(numCells);
  assert(meshDim);
  assert(_mesh);

#if 0 // Sieve
  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& meshCells = sieveMesh->heightStratum(0);
  assert(!meshCells.isNull());
  const SieveMesh::label_sequence::iterator meshCellsBegin = meshCells->begin();
  const SieveMesh::label_sequence::iterator meshCellsEnd = meshCells->end();

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
    const int coneSize = pV.getSize();
    assert(coneSize == *numCorners);
    for(int p = 0; p < coneSize; ++p, ++i) {
      (*cells)[i] = vNumbering->getIndex(cone[p]);
    }
    pV.clear();
  } // for
#endif // Sieve

  PetscDM dmMesh = _mesh->dmMesh();assert(dmMesh);
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  
  *numCells = _mesh->numCells();
  *numCorners = _mesh->coneSize();
  *meshDim = _mesh->dimension();
  assert(cellsStratum.size() == *numCells);

  cells->resize((*numCells)*(*numCorners));

  PetscIS globalVertexNumbers = NULL;
  const PetscInt* gvertex = NULL;
  const PetscInt* cone = NULL;
  PetscInt coneSize = 0, v = 0, count = 0;
  PetscErrorCode err = 0;

  err = DMPlexGetVertexNumbering(dmMesh, &globalVertexNumbers);PYLITH_CHECK_ERROR(err);
  err = ISGetIndices(globalVertexNumbers, &gvertex);PYLITH_CHECK_ERROR(err);
  for (PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    err = DMPlexGetConeSize(dmMesh, c, &coneSize);PYLITH_CHECK_ERROR(err);
    assert(coneSize == *numCorners);

    err = DMPlexGetCone(dmMesh, c, &cone);PYLITH_CHECK_ERROR(err);
    for(PetscInt p = 0; p < coneSize; ++p) {
      const PetscInt gv = gvertex[cone[p]-vStart];
      (*cells)[index++] = gv < 0 ? -(gv+1) : gv;
    } // for
  } // for  

  PYLITH_METHOD_END;
} // _getCells

// ----------------------------------------------------------------------
// Tag cells in mesh with material identifiers.
void
pylith::meshio::MeshIO::_setMaterials(const int_array& materialIds)
{ // _setMaterials
  PYLITH_METHOD_BEGIN;

  assert(_mesh);
  PetscErrorCode err;

  if (!_mesh->commRank()) {
#if 1 // Sieve
    const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
    if (!sieveMesh.isNull()) {
      const ALE::Obj<SieveMesh::label_type>& labelMaterials = sieveMesh->createLabel("material-id");
      const ALE::Obj<SieveMesh::label_sequence>& cells = sieveMesh->heightStratum(0);
      assert(!cells.isNull());
      const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
      const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

      const int numCells = materialIds.size();
      if (cells->size() != numCells) {
        std::ostringstream msg;
        msg << "Mismatch in size of materials identifier array ("
            << numCells << ") and number of cells in mesh ("<< cells->size() << ").";
        throw std::runtime_error(msg.str());
      } // if
      int i = 0;
      for(SieveMesh::label_sequence::iterator e_iter = cellsBegin; e_iter != cellsEnd; ++e_iter) {
        sieveMesh->setValue(labelMaterials, *e_iter, materialIds[i++]);
      }
    }
#endif // Sieve

    PetscDM dmMesh = _mesh->dmMesh();assert(dmMesh);
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();

    if (cellsStratum.size() != materialIds.size()) {
      std::ostringstream msg;
      msg << "Mismatch in size of materials identifier array ("
          << materialIds.size() << ") and number of cells in mesh ("<< (cEnd - cStart) << ").";
      throw std::runtime_error(msg.str());
    } // if
    PetscErrorCode err = 0;
    for(PetscInt c = cStart; c < cEnd; ++c) {
      err = DMPlexSetLabelValue(dmMesh, "material-id", c, materialIds[c-cStart]);PYLITH_CHECK_ERROR(err);
    } // for
  } // if

  PYLITH_METHOD_END;
} // _setMaterials

// ----------------------------------------------------------------------
// Get material identifiers for cells.
void
pylith::meshio::MeshIO::_getMaterials(int_array* materialIds) const
{ // _getMaterials
  PYLITH_METHOD_BEGIN;

  assert(materialIds);
  assert(_mesh);

#if 0 // Sieve
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
#endif // Sieve

  PetscDM dmMesh = _mesh->dmMesh();assert(dmMesh);
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  materialIds->resize(cellsStratum.size());
  PetscErrorCode err = 0;
  PetscInt matId = 0;
  for(PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    err = DMPlexGetLabelValue(dmMesh, "material-id", c, &matId);PYLITH_CHECK_ERROR(err);
    (*materialIds)[index++] = matId;
  } // for

  PYLITH_METHOD_END;
} // _getMaterials

// ----------------------------------------------------------------------
// Build a point group as an int section.
void
pylith::meshio::MeshIO::_setGroup(const std::string& name,
				  const GroupPtType type,
				  const int_array& points)
{ // _setGroup
  PYLITH_METHOD_BEGIN;

  assert(_mesh);

#if 1 // Sieve
  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  if (!sieveMesh.isNull()) {
    if (sieveMesh->hasIntSection(name)) {
      std::ostringstream msg;
      msg << "Could not setup group '" << name
          << "'. Group already exists in mesh.";
      throw std::runtime_error(msg.str());
    } // if

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
      groupField->setChart(IntSection::chart_type(numCells, numCells+numVertices));
      for(int i=0; i < numPoints; ++i)
        groupField->setFiberDimension(numCells+points[i], 1);
    } // if/else
    sieveMesh->allocate(groupField);
  }
#endif // Sieve

  PetscDM dmMesh = _mesh->dmMesh();assert(dmMesh);
  const PetscInt numPoints = points.size();
  PetscErrorCode err;

  if (CELL == type) {
    for(PetscInt p = 0; p < numPoints; ++p) {
      err = DMPlexSetLabelValue(dmMesh, name.c_str(), points[p], 1);PYLITH_CHECK_ERROR(err);
    } // for
  } else if (VERTEX == type) {
    PetscInt cStart, cEnd, numCells;

    err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    numCells = cEnd - cStart;
    for(PetscInt p = 0; p < numPoints; ++p) {
      err = DMPlexSetLabelValue(dmMesh, name.c_str(), numCells+points[p], 1);PYLITH_CHECK_ERROR(err);
    } // for
  } // if/else

  PYLITH_METHOD_END;
} // _setGroup

// ----------------------------------------------------------------------
// Create empty groups on other processes
void
pylith::meshio::MeshIO::_distributeGroups()
{ // _distributeGroups
  PYLITH_METHOD_BEGIN;

  assert(_mesh);

  // dmMesh does not need to broadcast the label names. They come from proc 0.

#if 1 // Sieve
  const ALE::Obj<SieveMesh>& sieveMesh = _mesh->sieveMesh();
  if (!sieveMesh.isNull()) {
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
  }
#endif // Sieve

  PYLITH_METHOD_END;
} // _distributeGroups

// ----------------------------------------------------------------------
// Get names of all groups in mesh.
void
pylith::meshio::MeshIO::_getGroupNames(string_vector* names) const
{ // _getGroups
  PYLITH_METHOD_BEGIN;

  assert(names);
  assert(_mesh);

#if 0 // Sieve
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
#endif // Sieve

  PetscDM dmMesh = _mesh->dmMesh();assert(dmMesh);
  PetscInt numGroups = 0;
  PetscErrorCode err = 0;
  err = DMPlexGetNumLabels(dmMesh, &numGroups);PYLITH_CHECK_ERROR(err);
  numGroups -= 2; // Remove depth and material labels.
  names->resize(numGroups);

  for (int iGroup=0, iLabel=numGroups-1; iGroup < numGroups; ++iGroup, --iLabel) {
    const char* namestr = NULL;
    err = DMPlexGetLabelName(dmMesh, iLabel, &namestr);PYLITH_CHECK_ERROR(err);
    (*names)[iGroup] = namestr;
  } // for

  PYLITH_METHOD_END;
} // _getGroups

// ----------------------------------------------------------------------
// Get group entities
void
pylith::meshio::MeshIO::_getGroup(int_array* points,
				  GroupPtType* groupType,
				  const char *name) const
{ // _getGroup
  PYLITH_METHOD_BEGIN;

  assert(points);
  assert(groupType);
  assert(_mesh);

#if 0
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
    *groupType = CELL;
    numbering = sieveMesh->getFactory()->getNumbering(sieveMesh, 
						      sieveMesh->depth());
  } else {
    *groupType = VERTEX;
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
#endif

  PetscDM dmMesh = _mesh->dmMesh();assert(dmMesh);
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  const PetscInt numCells = cellsStratum.size();

  PetscInt pStart, pEnd, firstPoint = 0;
  PetscErrorCode err = 0;
  err = DMPlexGetChart(dmMesh, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
  for(PetscInt p = pStart; p < pEnd; ++p) {
    PetscInt val;
    err = DMPlexGetLabelValue(dmMesh, name, p, &val);PYLITH_CHECK_ERROR(err);
    if (val >= 0) {
      firstPoint = p;
      break;
    } // if
  } // for
  *groupType = (firstPoint >= cStart && firstPoint < cEnd) ? CELL : VERTEX;

  PetscInt groupSize;
  err = DMPlexGetStratumSize(dmMesh, name, 1, &groupSize);PYLITH_CHECK_ERROR(err);
  points->resize(groupSize);

  const PetscInt offset = (VERTEX == *groupType) ? numCells : 0;
  PetscIS groupIS = NULL;
  const PetscInt* groupIndices = NULL;
  err = DMPlexGetStratumIS(dmMesh, name, 1, &groupIS);PYLITH_CHECK_ERROR(err);
  err = ISGetIndices(groupIS, &groupIndices);PYLITH_CHECK_ERROR(err);
  for(PetscInt p = 0; p < groupSize; ++p) {
    (*points)[p] = groupIndices[p]-offset;
  } // for
  err = ISRestoreIndices(groupIS, &groupIndices);PYLITH_CHECK_ERROR(err);
  err = ISDestroy(&groupIS);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _getGroup


// End of file
