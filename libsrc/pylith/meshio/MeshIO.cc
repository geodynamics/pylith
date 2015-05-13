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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "MeshIO.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "pylith/utils/array.hh" // USES scalar_array, int_array
#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

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

  // Check mesh consistency
  topology::MeshOps::checkTopology(*_mesh);
  // Respond to PETSc diagnostic output
  PetscErrorCode err = DMViewFromOptions(_mesh->dmMesh(), "pylith_", "-dm_view");PYLITH_CHECK_ERROR(err);

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

  PetscDM dmMesh = _mesh->dmMesh();assert(dmMesh);
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  
  *numCells = _mesh->numCells();
  *numCorners = _mesh->numCorners();
  *meshDim = _mesh->dimension();
  assert(cellsStratum.size() == *numCells);

  cells->resize((*numCells)*(*numCorners));

  PetscIS globalVertexNumbers = NULL;
  const PetscInt* gvertex = NULL;
  PetscErrorCode err = 0;

  err = DMPlexGetVertexNumbering(dmMesh, &globalVertexNumbers);PYLITH_CHECK_ERROR(err);
  err = ISGetIndices(globalVertexNumbers, &gvertex);PYLITH_CHECK_ERROR(err);
  for (PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    PetscInt nC = 0, closureSize, *closure = NULL;

    err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    for (PetscInt cl = 0; cl < closureSize*2; cl += 2) {
      if ((closure[cl] >= vStart) && (closure[cl] < vEnd)) {
        const PetscInt gv = gvertex[closure[cl]-vStart];
        (*cells)[index++] = gv < 0 ? -(gv+1) : gv;
        ++nC;
      }
    } // for
    err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    err = DMPlexInvertCell(*meshDim, nC, &(*cells)[index-nC]);PYLITH_CHECK_ERROR(err);
    assert(nC == *numCorners);
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

  if (!_mesh->commRank()) {
    PetscDM dmMesh = _mesh->dmMesh();assert(dmMesh);
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();

    if (size_t(cellsStratum.size()) != materialIds.size()) {
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

  PetscDM        dmMesh    = _mesh->dmMesh();assert(dmMesh);
  const PetscInt numPoints = points.size();
  DMLabel        label;
  PetscErrorCode err;

  err = DMPlexCreateLabel(dmMesh, name.c_str());PYLITH_CHECK_ERROR(err);
  err = DMPlexGetLabel(dmMesh, name.c_str(), &label);PYLITH_CHECK_ERROR(err);
  if (CELL == type) {
    for(PetscInt p = 0; p < numPoints; ++p) {
      err = DMLabelSetValue(label, points[p], 1);PYLITH_CHECK_ERROR(err);
    } // for
  } else if (VERTEX == type) {
    PetscInt cStart, cEnd, vStart, vEnd, numCells;

    err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
    numCells = cEnd - cStart;
    for(PetscInt p = 0; p < numPoints; ++p) {
      err = DMLabelSetValue(label, numCells+points[p], 1);PYLITH_CHECK_ERROR(err);
    } // for
    // Also add any non-cells which have all vertices marked
    for(PetscInt p = 0; p < numPoints; ++p) {
      const PetscInt vertex = numCells+points[p];
      PetscInt      *star   = NULL, starSize, s;

      err = DMPlexGetTransitiveClosure(dmMesh, vertex, PETSC_FALSE, &starSize, &star);PYLITH_CHECK_ERROR(err);
      for (s = 0; s < starSize*2; s += 2) {
        const PetscInt point   = star[s];
        PetscInt      *closure = NULL, closureSize, c, value;
        PetscBool      marked  = PETSC_TRUE;

        if ((point >= cStart) && (point < cEnd)) continue;
        err = DMPlexGetTransitiveClosure(dmMesh, point, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        for (c = 0; c < closureSize*2; c += 2) {
          if ((closure[c] >= vStart) && (closure[c] < vEnd)) {
            err = DMLabelGetValue(label, closure[c], &value);PYLITH_CHECK_ERROR(err);
            if (value != 1) {marked = PETSC_FALSE; break;}
          }
        }
        err = DMPlexRestoreTransitiveClosure(dmMesh, point, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        if (marked) {err = DMLabelSetValue(label, point, 1);PYLITH_CHECK_ERROR(err);}
      }
      err = DMPlexRestoreTransitiveClosure(dmMesh, vertex, PETSC_FALSE, &starSize, &star);PYLITH_CHECK_ERROR(err);
    }
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

  PetscDM dmMesh = _mesh->dmMesh();assert(dmMesh);
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  const PetscInt numCells = cellsStratum.size();

  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  PetscIS groupIS = NULL;
  const PetscInt* groupIndices = NULL;
  PetscErrorCode err;
  err = DMPlexGetStratumIS(dmMesh, name, 1, &groupIS);PYLITH_CHECK_ERROR(err);
  err = ISGetIndices(groupIS, &groupIndices);PYLITH_CHECK_ERROR(err);

  PetscInt totalSize;
  err = DMPlexGetStratumSize(dmMesh, name, 1, &totalSize);PYLITH_CHECK_ERROR(err);

  *groupType = VERTEX;
  if (totalSize > 0 && (groupIndices[0] >= cStart && groupIndices[0] < cEnd)) {
    *groupType = CELL;
  } // if
    
  PetscInt offset = 0;
  PetscInt pStart = cStart;
  PetscInt pEnd = cEnd;
  if (VERTEX == *groupType) {
    offset = numCells;
    pStart = vStart;
    pEnd = vEnd;
  } // if

  // Count number of cells/vertices, filtering out edges and faces
  PetscInt groupSize = 0;
  for (PetscInt p = 0; p < totalSize; ++p) {
    if (groupIndices[p] >= pStart && groupIndices[p] < pEnd) {
      ++groupSize;
    } // if
  } // for

  points->resize(groupSize);
  for(PetscInt p = 0; p < groupSize; ++p) {
    if (groupIndices[p] >= pStart && groupIndices[p] < pEnd) {
      (*points)[p] = groupIndices[p]-offset;
    } // if
  } // for
  err = ISRestoreIndices(groupIS, &groupIndices);PYLITH_CHECK_ERROR(err);
  err = ISDestroy(&groupIS);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _getGroup


// End of file
