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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Field.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
template<typename mesh_type>
pylith::topology::Field<mesh_type>::Field(const mesh_type& mesh) :
  _mesh(mesh)
{ // constructor
  _metadata["default"].label = "unknown";
  _metadata["default"].vectorFieldType = OTHER;
  _metadata["default"].scale = 1.0;
  _metadata["default"].dimsOkay = false;
  if (mesh.dmMesh()) {
    PetscDM dm = mesh.dmMesh();
    PetscVec coordVec = NULL;
    PetscSection s = NULL;
    PetscErrorCode err;

    err = DMPlexClone(dm, &_dm);CHECK_PETSC_ERROR(err);
    err = DMGetCoordinatesLocal(dm, &coordVec);CHECK_PETSC_ERROR(err);
    if (coordVec) {
      PetscDM coordDM=NULL, newCoordDM=NULL;
      PetscSection coordSection=NULL, newCoordSection=NULL;

      err = DMGetCoordinateDM(dm, &coordDM);CHECK_PETSC_ERROR(err);
      err = DMGetCoordinateDM(_dm, &newCoordDM);CHECK_PETSC_ERROR(err);
      err = DMGetDefaultSection(coordDM, &coordSection);CHECK_PETSC_ERROR(err);
      err = PetscSectionClone(coordSection, &newCoordSection);CHECK_PETSC_ERROR(err);
      err = DMSetDefaultSection(newCoordDM, newCoordSection);CHECK_PETSC_ERROR(err);
      err = DMSetCoordinatesLocal(_dm, coordVec);CHECK_PETSC_ERROR(err);
    }
    err = PetscSectionCreate(mesh.comm(), &s);CHECK_PETSC_ERROR(err);
    err = DMSetDefaultSection(_dm, s);CHECK_PETSC_ERROR(err);
  } else {
    _dm = NULL;
  }
  _globalVec = NULL;
  _localVec  = NULL;
} // constructor

// ----------------------------------------------------------------------
// Constructor with mesh, DM, and metadata
template<typename mesh_type>
pylith::topology::Field<mesh_type>::Field(const mesh_type& mesh,
					  PetscDM dm,
					  const Metadata& metadata) :
  _mesh(mesh),
  _dm(dm)
{ // constructor
  assert(dm);
  PetscErrorCode err;

  _metadata["default"] = metadata;
  err = DMCreateLocalVector(_dm, &_localVec);CHECK_PETSC_ERROR(err);
  err = DMCreateGlobalVector(_dm, &_globalVec);CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) _globalVec, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) _localVec,  _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
} // constructor

// ----------------------------------------------------------------------
// Constructor with mesh, DM, local data, and metadata
template<typename mesh_type>
pylith::topology::Field<mesh_type>::Field(const mesh_type& mesh,
					  PetscDM dm,
					  PetscVec localVec,
					  const Metadata& metadata) :
  _mesh(mesh),
  _dm(dm)
{ // constructor
  assert(dm);
  assert(localVec);

  PetscErrorCode err;

  _metadata["default"] = metadata;
  err = DMCreateLocalVector(_dm, &_localVec);CHECK_PETSC_ERROR(err);
  err = DMCreateGlobalVector(_dm, &_globalVec);CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) _globalVec, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) _localVec,  _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  err = VecCopy(localVec, _localVec);CHECK_PETSC_ERROR(err);
} // constructor

// ----------------------------------------------------------------------
// Constructor with field and subfields
template<typename mesh_type>
pylith::topology::Field<mesh_type>::Field(const Field& src,
					  const int fields[],
					  int numFields) :
  _mesh(src._mesh)
{ // constructor
  PetscDM dm = mesh.dmMesh(), coordDM=NULL, newCoordDM=NULL;
  PetscSection coordSection=NULL, newCoordSection=NULL;
  PetscVec coordVec=NULL;
  PetscSection s=NULL;
  PetscErrorCode err;

  assert(dm);
  assert(src._dm);

  _metadata["default"] = src._metadata["default"];
  err = DMGetDefaultSection(src._dm, &s);CHECK_PETSC_ERROR(err);
  for(PetscInt f = 0; f < numFields; ++f) {
    const char *name;

    err = PetscSectionGetFieldName(s, fields[f], &name);CHECK_PETSC_ERROR(err);
    _metadata[name] = src._metadata[name];
  } // for
  err = DMCreateSubDM(dm, numFields, fields, NULL, &_dm);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dm, &coordVec);CHECK_PETSC_ERROR(err);
  if (coordVec) {
    err = DMGetCoordinateDM(dm, &coordDM);CHECK_PETSC_ERROR(err);
    err = DMGetCoordinateDM(_dm, &newCoordDM);CHECK_PETSC_ERROR(err);
    err = DMGetDefaultSection(coordDM, &coordSection);CHECK_PETSC_ERROR(err);
    err = PetscSectionClone(coordSection, &newCoordSection);CHECK_PETSC_ERROR(err);
    err = DMSetDefaultSection(newCoordDM, newCoordSection);CHECK_PETSC_ERROR(err);
    err = DMSetCoordinatesLocal(_dm, coordVec);CHECK_PETSC_ERROR(err);
  } // if
  _globalVec = NULL;
  _localVec  = NULL;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
template<typename mesh_type>
pylith::topology::Field<mesh_type>::~Field(void)
{ // destructor
  deallocate();
  PetscErrorCode err = DMDestroy(&_dm);CHECK_PETSC_ERROR(err);
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::deallocate(void)
{ // deallocate
  PetscErrorCode err = 0;
  
  const typename scatter_map_type::const_iterator scattersEnd = _scatters.end();
  for (typename scatter_map_type::iterator s_iter=_scatters.begin();
       s_iter != scattersEnd;
       ++s_iter) {

    err = DMDestroy(&s_iter->second.dm);CHECK_PETSC_ERROR(err);
    err = VecDestroy(&s_iter->second.vector);CHECK_PETSC_ERROR(err);

    if (s_iter->second.scatter) {
      err = VecScatterDestroy(&s_iter->second.scatter);CHECK_PETSC_ERROR(err);
    } // if
    err = VecDestroy(&s_iter->second.scatterVec);CHECK_PETSC_ERROR(err);
  } // for
  _scatters.clear();
  err = VecDestroy(&_globalVec);CHECK_PETSC_ERROR(err);
  err = VecDestroy(&_localVec);CHECK_PETSC_ERROR(err);
} // deallocate

// ----------------------------------------------------------------------
// Set label for field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::label(const char* value)
{ // label
  PetscErrorCode err;

  _metadata["default"].label = value;
  if (_localVec)  {err = PetscObjectSetName((PetscObject) _localVec, value);CHECK_PETSC_ERROR(err);}
  if (_globalVec) {err = PetscObjectSetName((PetscObject) _globalVec, value);CHECK_PETSC_ERROR(err);}

  const typename scatter_map_type::const_iterator scattersEnd = _scatters.end();
  for (typename scatter_map_type::const_iterator s_iter=_scatters.begin();
       s_iter != scattersEnd;
       ++s_iter) {
    if (s_iter->second.vector) {
      err = PetscObjectSetName((PetscObject)s_iter->second.vector, value);CHECK_PETSC_ERROR(err);    
    } // if
  } // for
} // label

// ----------------------------------------------------------------------
// Get spatial dimension of domain.
template<typename mesh_type>
int
pylith::topology::Field<mesh_type>::spaceDim(void) const
{ // spaceDim
  const spatialdata::geocoords::CoordSys* cs = _mesh.coordsys();
  return (cs) ? cs->spaceDim() : 0;
} // spaceDim

// ----------------------------------------------------------------------
// Get the chart size.
template<typename mesh_type>
int
pylith::topology::Field<mesh_type>::chartSize(void) const
{ // chartSize
  assert(_dm);
  PetscSection s = NULL;
  PetscInt pStart, pEnd;
  PetscErrorCode err;

  err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetChart(s, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
  return pEnd-pStart;
} // chartSize

// ----------------------------------------------------------------------
// Get the number of degrees of freedom.
template<typename mesh_type>
int
pylith::topology::Field<mesh_type>::sectionSize(void) const
{ // sectionSize
  PetscInt size = 0;

  if (_dm) {
    PetscSection s = NULL;
    PetscErrorCode err;

    err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetStorageSize(s, &size);CHECK_PETSC_ERROR(err);
  }
  return size;
} // sectionSize

// ----------------------------------------------------------------------
// Create seive section.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(void)
{ // newSection
  // Clear memory
  clear();
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion for a list of
// points.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(const int_array& points,
                                               const int fiberDim)
{ // newSection
  typedef typename mesh_type::SieveMesh::point_type point_type;
  PetscErrorCode err;

  // Clear memory
  clear();
  assert(_dm);
  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _metadata["default"].label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if
  
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");
  const PetscInt npts = points.size();
  if (npts > 0) {
    PetscSection s = NULL;
    PetscInt pointMin = 0, pointMax = 0;

    for (PetscInt i = 0; i < npts; ++i) {
      pointMin = std::min(pointMin, points[i]);
      pointMax = std::max(pointMax, points[i]);
    } // for
    err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
    err = DMSetDefaultGlobalSection(_dm, NULL);CHECK_PETSC_ERROR(err);
    err = PetscSectionSetChart(s, pointMin, pointMax+1);CHECK_PETSC_ERROR(err);
    for (PetscInt i = 0; i < npts; ++i) {
      err = PetscSectionSetDof(s, points[i], fiberDim);CHECK_PETSC_ERROR(err);
    } // for
  } else { // create empty chart
    PetscSection s = NULL;

    err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
    err = DMSetDefaultGlobalSection(_dm, NULL);CHECK_PETSC_ERROR(err);
    err = PetscSectionSetChart(s, 0, 0);CHECK_PETSC_ERROR(err);
  } // if/else

  logger.stagePop();
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion for a list of
// points.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(const PetscInt *points, const PetscInt num,
                                               const int fiberDim)
{ // newSection
  PetscErrorCode err;

  // Clear memory
  clear();
  assert(_dm);
  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _metadata["default"].label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if
  
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");
  if (num > 0) {
    PetscSection s = NULL;
    PetscInt pointMin = 0, pointMax = 0;

    for (PetscInt i = 0; i < num; ++i) {
      pointMin = std::min(pointMin, points[i]);
      pointMax = std::max(pointMax, points[i]);
    } // for
    err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
    err = DMSetDefaultGlobalSection(_dm, NULL);CHECK_PETSC_ERROR(err);
    err = PetscSectionSetChart(s, pointMin, pointMax+1);CHECK_PETSC_ERROR(err);
    for (PetscInt i = 0; i < num; ++i) {
      err = PetscSectionSetDof(s, points[i], fiberDim);CHECK_PETSC_ERROR(err);
    } // for
  } else { // create empty chart
    PetscSection s = NULL;

    err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
    err = DMSetDefaultGlobalSection(_dm, NULL);CHECK_PETSC_ERROR(err);
    err = PetscSectionSetChart(s, 0, 0);CHECK_PETSC_ERROR(err);
  } // if/else

  logger.stagePop();
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(const DomainEnum domain,
                                               const int fiberDim,
                                               const int stratum)
{ // newSection
  // Changing this because cells/vertices are numbered differently in the new scheme
  assert(_dm);
  PetscSection s = NULL;
  PetscInt pStart, pEnd;
  PetscErrorCode err;

  switch(domain) {
  case VERTICES_FIELD:
    err = DMPlexGetDepthStratum(_dm, stratum, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    break;
  case CELLS_FIELD:
    err = DMPlexGetHeightStratum(_dm, stratum, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    break;
  case FACES_FIELD:
    err = DMPlexGetHeightStratum(_dm, stratum+1, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    break;
  case POINTS_FIELD:
    err = DMPlexGetChart(_dm, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    break;
  default:
    std::cerr << "Unknown value for DomainEnum: " << domain << std::endl;
    throw std::logic_error("Bad domain enum in Field.");
  }
  newSection(pStart, pEnd, fiberDim);
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(const PetscInt pStart, const PetscInt pEnd,
                                               const int fiberDim)
{ // newSection
  // Changing this because cells/vertices are numbered differently in the new scheme
  assert(_dm);
  PetscSection s = NULL;
  PetscErrorCode err;

  err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
  err = DMSetDefaultGlobalSection(_dm, NULL);CHECK_PETSC_ERROR(err);
  err = PetscSectionSetChart(s, pStart, pEnd);CHECK_PETSC_ERROR(err);

  for(PetscInt p = pStart; p < pEnd; ++p) {
    err = PetscSectionSetDof(s, p, fiberDim);CHECK_PETSC_ERROR(err);
  }
} // newSection

// ----------------------------------------------------------------------
// Create section given chart.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(const Field& src,
                                               const int fiberDim)
{ // newSection
  // Clear memory
  clear();
  assert(_dm);assert(src._dm);

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");

  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _metadata["default"].label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  PetscSection srcs=NULL, s=NULL;
  PetscInt pStart, pEnd;
  PetscErrorCode err;

  err = DMGetDefaultSection(src._dm, &srcs);CHECK_PETSC_ERROR(err);
  err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
  err = DMSetDefaultGlobalSection(_dm, NULL);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetChart(srcs, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
  err = PetscSectionSetChart(s, pStart, pEnd);CHECK_PETSC_ERROR(err);
  for(PetscInt p = pStart; p < pEnd; ++p) {
    err = PetscSectionSetDof(s, p, fiberDim);CHECK_PETSC_ERROR(err);
  }

  logger.stagePop();
} // newSection

// ----------------------------------------------------------------------
// Create section with same layout as another section.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::cloneSection(const Field& src)
{ // cloneSection
  std::string origLabel = _metadata["default"].label;

  // Clear memory
  clear();

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");

  _metadata["default"] = const_cast<Field&>(src)._metadata["default"];
  label(origLabel.c_str());

  PetscSection section = src.petscSection();
  PetscSection newSection = NULL;
  PetscErrorCode err;

  assert(_dm);
  err = PetscSectionClone(section, &newSection);CHECK_PETSC_ERROR(err);
  err = DMSetDefaultSection(_dm, newSection);CHECK_PETSC_ERROR(err);
  err = DMCreateGlobalVector(_dm, &_globalVec);CHECK_PETSC_ERROR(err);
  err = DMCreateLocalVector(_dm, &_localVec);CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) _globalVec, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) _localVec,  _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
    
  // Reuse scatters in clone
  if (src._scatters.size() > 0) {
    const typename scatter_map_type::const_iterator scattersEnd = src._scatters.end();
    for (typename scatter_map_type::const_iterator s_iter=src._scatters.begin();
         s_iter != scattersEnd;
         ++s_iter) {
      ScatterInfo sinfo;
      sinfo.vector = 0;
      sinfo.scatter = 0;
      sinfo.scatterVec = 0;

      // Copy DM
      sinfo.dm = s_iter->second.dm;
      err = PetscObjectReference((PetscObject) sinfo.dm);
      CHECK_PETSC_ERROR(err);

      // Copy scatter
      if (s_iter->second.scatter) {
        sinfo.scatter = s_iter->second.scatter;
        err = PetscObjectReference((PetscObject) sinfo.scatter);
        CHECK_PETSC_ERROR(err);
      } // if
      
      // Create scatter Vec
      sinfo.scatterVec = _localVec;
      err = PetscObjectReference((PetscObject) sinfo.scatterVec);CHECK_PETSC_ERROR(err);

      // Create vector using sizes from source section
      PetscInt vecLocalSize = 0;
      PetscInt vecGlobalSize = 0, vecGlobalSize2 = 0;
      err = VecGetLocalSize(s_iter->second.vector, &vecLocalSize);CHECK_PETSC_ERROR(err);
      err = VecGetSize(s_iter->second.vector, &vecGlobalSize);CHECK_PETSC_ERROR(err);
      err = VecGetSize(_globalVec, &vecGlobalSize2);CHECK_PETSC_ERROR(err);
      
      if (vecGlobalSize != vecGlobalSize2) {
        int blockSize = 1;
        err = VecGetBlockSize(s_iter->second.vector, &blockSize);CHECK_PETSC_ERROR(err);
        err = VecCreate(_mesh.comm(), &sinfo.vector);CHECK_PETSC_ERROR(err);
        err = VecSetSizes(sinfo.vector, vecLocalSize, vecGlobalSize);CHECK_PETSC_ERROR(err);
        err = VecSetBlockSize(sinfo.vector, blockSize);CHECK_PETSC_ERROR(err);
        err = VecSetFromOptions(sinfo.vector); CHECK_PETSC_ERROR(err);  
      } else {
        sinfo.vector = _globalVec;
        err = PetscObjectReference((PetscObject) sinfo.vector);
        CHECK_PETSC_ERROR(err);
      } // if/else
      err = PetscObjectSetName((PetscObject)sinfo.vector, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
	
      _scatters[s_iter->first] = sinfo;
    } // for
  } // if
  logger.stagePop();
} // cloneSection

// ----------------------------------------------------------------------
// Clear variables associated with section.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::clear(void)
{ // clear
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");

  deallocate();
  _metadata["default"].scale = 1.0;
  _metadata["default"].vectorFieldType = OTHER;
  _metadata["default"].dimsOkay = false;

  logger.stagePop();
} // clear

// ----------------------------------------------------------------------
// Allocate Sieve section.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::allocate(void)
{ // allocate
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");
  PetscSection   s = NULL;
  PetscErrorCode err;

  if (_dm) {
    err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
  } // if
  assert(s);

  err = PetscSectionSetUp(s);CHECK_PETSC_ERROR(err);
  err = DMCreateGlobalVector(_dm, &_globalVec);CHECK_PETSC_ERROR(err);
  err = DMCreateLocalVector(_dm, &_localVec);CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) _globalVec, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) _localVec,  _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);

  logger.stagePop();
} // allocate

// ----------------------------------------------------------------------
// Zero section values (excluding constrained DOF).
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::zero(void)
{ // zero
  assert(_localVec);
  PetscSection section = NULL;
  PetscInt pStart, pEnd, maxDof = 0;
  PetscErrorCode err;

  err = DMGetDefaultSection(_dm, &section);CHECK_PETSC_ERROR(err);    
  err = PetscSectionGetChart(section, &pStart, &pEnd);CHECK_PETSC_ERROR(err);    
  if (pEnd > pStart) {err = PetscSectionGetMaxDof(section, &maxDof);CHECK_PETSC_ERROR(err);}
  scalar_array values(maxDof);
  values *= 0.0;

  for(PetscInt p = pStart; p < pEnd; ++p) {
    PetscInt dof;

    err = PetscSectionGetDof(section, p, &dof);CHECK_PETSC_ERROR(err);
    if (dof > 0) {
      assert(dof <= maxDof);
      err = DMPlexVecSetClosure(_dm, section, _localVec, p, &values[0], INSERT_VALUES);CHECK_PETSC_ERROR(err);
    } // if
  } // for
} // zero

// ----------------------------------------------------------------------
// Zero section values (including constrained DOF).
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::zeroAll(void)
{ // zeroAll
  assert(_localVec);
  PetscErrorCode err = VecSet(_localVec, 0.0);CHECK_PETSC_ERROR(err);
} // zeroAll

// ----------------------------------------------------------------------
// Complete section by assembling across processors.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::complete(void)
{ // complete
  assert(_dm);
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Completion");
  // Not sure if DMLocalToLocal() would work
  PetscErrorCode err;

  err = VecSet(_globalVec, 0.0);CHECK_PETSC_ERROR(err);
  err = DMLocalToGlobalBegin(_dm, _localVec, ADD_VALUES, _globalVec);CHECK_PETSC_ERROR(err);
  err = DMLocalToGlobalEnd(_dm, _localVec, ADD_VALUES, _globalVec);CHECK_PETSC_ERROR(err);
  err = DMGlobalToLocalBegin(_dm, _globalVec, INSERT_VALUES, _localVec);CHECK_PETSC_ERROR(err);
  err = DMGlobalToLocalEnd(_dm, _globalVec, INSERT_VALUES, _localVec);CHECK_PETSC_ERROR(err);
  logger.stagePop();
} // complete

// ----------------------------------------------------------------------
// Copy field values and metadata.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::copy(const Field& field)
{ // copy
  // Check compatibility of sections
  const int srcSize = field.chartSize();
  const int dstSize = chartSize();
  if (field.spaceDim() != spaceDim() ||
      const_cast<Field&>(field)._metadata["default"].vectorFieldType != _metadata["default"].vectorFieldType ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot copy values from section '" << const_cast<Field&>(field)._metadata["default"].label 
	<< "' to section '" << _metadata["default"].label
	<< "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    space dim: " << field.spaceDim() << "\n"
	<< "    vector field type: " << const_cast<Field&>(field)._metadata["default"].vectorFieldType << "\n"
	<< "    scale: " << const_cast<Field&>(field)._metadata["default"].scale << "\n"
	<< "    size: " << srcSize << "\n"
	<< "  Destination section:\n"
	<< "    space dim: " << spaceDim() << "\n"
	<< "    vector field type: " << _metadata["default"].vectorFieldType << "\n"
	<< "    scale: " << _metadata["default"].scale << "\n"
	<< "    size: " << dstSize;
    throw std::runtime_error(msg.str());
  } // if
  assert(_localVec && field._localVec);

  PetscErrorCode err = VecCopy(field._localVec, _localVec);CHECK_PETSC_ERROR(err);

  label(const_cast<Field&>(field)._metadata["default"].label.c_str()); // Update label
  _metadata["default"].scale = const_cast<Field&>(field)._metadata["default"].scale;
} // copy

template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::copy(PetscSection osection,
					 PetscInt field,
					 PetscInt component,
					 PetscVec ovec)
{ // copy
  assert(osection);
  assert(ovec);
  assert(_localVec);

  PetscSection section = NULL;
  PetscScalar *array = NULL, *oarray = NULL;
  PetscInt numFields, numComp, pStart, pEnd, qStart, qEnd;
  PetscErrorCode err;

  assert(_dm);
  err = DMGetDefaultSection(_dm, &section);CHECK_PETSC_ERROR(err);assert(section);  
  err = PetscSectionGetNumFields(osection, &numFields);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetChart(section, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetChart(osection, &qStart, &qEnd);CHECK_PETSC_ERROR(err);
  if ((pStart != qStart) || (pEnd != qEnd)) {
    std::ostringstream msg;

    msg << "Cannot copy values from Sieve section "
	<< _metadata["default"].label << "'. Sections are incompatible.\n"
	<< "  Source section:\n"
    << "    chart: ["<<pStart<<","<<pEnd<<")\n"
	<< "  Destination section:\n"
	<< "    space dim: " << spaceDim() << "\n"
	<< "    vector field type: " << _metadata["default"].vectorFieldType << "\n"
	<< "    scale: " << _metadata["default"].scale << "\n"
    << "    chart: ["<<qStart<<","<<qEnd<<")";
    throw std::runtime_error(msg.str());
  } // if
  if (field >= numFields) {
    std::ostringstream msg;
    msg << "Invalid field "<<field<<" should be in [0, "<<numFields<<")";
    throw std::runtime_error(msg.str());
  }
  if ((field >= 0) && (component >= 0)) {
    err = PetscSectionGetFieldComponents(osection, field, &numComp);CHECK_PETSC_ERROR(err);
    if (component >= numComp) {
      std::ostringstream msg;
      msg << "Invalid field component "<<component<<" should be in [0, "<<numComp<<")";
      throw std::runtime_error(msg.str());
    } // if
  } // if
  // Copy values from field
  err = VecGetArray(_localVec, &array);CHECK_PETSC_ERROR(err);
  err = VecGetArray(ovec, &oarray);CHECK_PETSC_ERROR(err);
  for(PetscInt p = pStart; p < pEnd; ++p) {
    PetscInt dof, off, odof, ooff;

    err = PetscSectionGetDof(section, p, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(section, p, &off);CHECK_PETSC_ERROR(err);
    if (field >= 0) {
      err = PetscSectionGetFieldDof(osection, p, field, &odof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetFieldOffset(osection, p, field, &ooff);CHECK_PETSC_ERROR(err);
      if (component >= 0) {
        assert(!(odof%numComp));
        odof  = odof/numComp;
        ooff += odof*component;
      } // if
    } else {
      err = PetscSectionGetDof(osection, p, &odof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(osection, p, &ooff);CHECK_PETSC_ERROR(err);
    } // else
    assert(odof == dof);
    if (!odof) continue;
    for(PetscInt d = 0; d < dof; ++d) {
      array[off+d] = oarray[ooff+d];
    } // for
  } // for
  err = VecRestoreArray(_localVec, &array);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(ovec, &oarray);CHECK_PETSC_ERROR(err);
} // copy

// ----------------------------------------------------------------------
// Add two fields, storing the result in one of the fields.
template<typename mesh_type>
pylith::topology::Field<mesh_type>&
pylith::topology::Field<mesh_type>::operator+=(const Field& field)
{ // operator+=
  // Check compatibility of sections
  const int srcSize = field.chartSize();
  const int dstSize = chartSize();
  if (field.spaceDim() != spaceDim() ||
      const_cast<Field&>(field)._metadata["default"].vectorFieldType != _metadata["default"].vectorFieldType ||
      const_cast<Field&>(field)._metadata["default"].scale != _metadata["default"].scale ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot add values from section '" << const_cast<Field&>(field)._metadata["default"].label 
	<< "' to section '" << _metadata["default"].label
	<< "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    space dim: " << field.spaceDim() << "\n"
    << "    vector field type: " << const_cast<Field&>(field)._metadata["default"].vectorFieldType << "\n"
    << "    scale: " << const_cast<Field&>(field)._metadata["default"].scale << "\n"
	<< "    size: " << srcSize << "\n"
	<< "  Destination section:\n"
	<< "    space dim: " << spaceDim() << "\n"
	<< "    vector field type: " << _metadata["default"].vectorFieldType << "\n"
	<< "    scale: " << _metadata["default"].scale << "\n"
	<< "    size: " << dstSize;
    throw std::runtime_error(msg.str());
  } // if
  assert(_localVec && field._localVec);
  PetscErrorCode err = VecAXPY(_localVec, 1.0, field._localVec);CHECK_PETSC_ERROR(err);
  return *this;
} // operator+=

// ----------------------------------------------------------------------
// Dimensionalize field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::dimensionalize(void) const
{ // dimensionalize
  if (!const_cast<Field*>(this)->_metadata["default"].dimsOkay) {
    std::ostringstream msg;
    msg << "Cannot dimensionalize field '" << const_cast<Field*>(this)->_metadata["default"].label
	<< "' because the flag "
	<< "has been set to keep field nondimensional.";
    throw std::runtime_error(msg.str());
  } // if

  spatialdata::units::Nondimensional normalizer;
  assert(_localVec);
  PetscSection section = NULL;
  PetscScalar *array = NULL;
  PetscInt pStart, pEnd;
  PetscErrorCode err;

  err = DMGetDefaultSection(_dm, &section);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetChart(section, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
  err = VecGetArray(_localVec, &array);CHECK_PETSC_ERROR(err);
  for(PetscInt p = pStart; p < pEnd; ++p) {
    PetscInt dof, off;

    err = PetscSectionGetDof(section, p, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(section, p, &off);CHECK_PETSC_ERROR(err);
    if (dof) {
      normalizer.dimensionalize(&array[off], dof, const_cast<Field*>(this)->_metadata["default"].scale);
    }
  }
  err = VecRestoreArray(_localVec, &array);CHECK_PETSC_ERROR(err);
} // dimensionalize

// ----------------------------------------------------------------------
// Print field to standard out.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::view(const char* label) const
{ // view
  std::string vecFieldString;
  switch(const_cast<Field*>(this)->_metadata["default"].vectorFieldType)
    { // switch
    case SCALAR:
      vecFieldString = "scalar";
      break;
    case VECTOR:
      vecFieldString = "vector";
      break;
    case TENSOR:
      vecFieldString = "tensor";
      break;
    case OTHER:
      vecFieldString = "other";
      break;
    case MULTI_SCALAR:
      vecFieldString = "multiple scalars";
      break;
    case MULTI_VECTOR:
      vecFieldString = "multiple vectors";
      break;
    case MULTI_TENSOR:
      vecFieldString = "multiple tensors";
      break;
    case MULTI_OTHER:
      vecFieldString = "multiple other values";
      break;
    default :
      std::cerr << "Unknown vector field value '" << const_cast<Field*>(this)->_metadata["default"].vectorFieldType
		<< "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad vector field type in Field.");
    } // switch

  std::cout << "Viewing field '" << const_cast<Field*>(this)->_metadata["default"].label << "' "<< label << ".\n"
	    << "  vector field type: " << vecFieldString << "\n"
	    << "  scale: " << const_cast<Field*>(this)->_metadata["default"].scale << "\n"
	    << "  dimensionalize flag: " << const_cast<Field*>(this)->_metadata["default"].dimsOkay << std::endl;
  if (_dm) {
    PetscSection section = NULL;
    PetscErrorCode err;

    err = DMGetDefaultSection(_dm, &section);CHECK_PETSC_ERROR(err);
    err = DMView(_dm, PETSC_VIEWER_STDOUT_WORLD);CHECK_PETSC_ERROR(err);
    err = PetscSectionView(section, PETSC_VIEWER_STDOUT_WORLD);CHECK_PETSC_ERROR(err);
    err = VecView(_localVec, PETSC_VIEWER_STDOUT_WORLD);CHECK_PETSC_ERROR(err);
  }
} // view

// ----------------------------------------------------------------------
// Create PETSc vector scatter for field. This is used to transfer
// information from the "global" PETSc vector view to the "local"
// Sieve section view.
template<typename mesh_type>
template<typename scatter_mesh_type>
void
pylith::topology::Field<mesh_type>::createScatter(const scatter_mesh_type& mesh,
						  const char* context)
{ // createScatter
  assert(context);
  PetscErrorCode err = 0;

  const bool createScatterOk = true;
  ScatterInfo& sinfo = _getScatter(context, createScatterOk);
  if (sinfo.scatter) {
    assert(sinfo.scatterVec);
    assert(sinfo.vector);
    return;
  } // if

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("GlobalOrder");
  PetscInt localSize, globalSize;

  err = PetscObjectReference((PetscObject) _dm);CHECK_PETSC_ERROR(err);
  err = PetscObjectReference((PetscObject) _globalVec);CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) _globalVec, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  err = VecGetSize(_localVec,  &localSize);CHECK_PETSC_ERROR(err);
  err = VecGetSize(_globalVec, &globalSize);CHECK_PETSC_ERROR(err);
  //assert(order->getLocalSize()  == localSize);
  //assert(order->getGlobalSize() == globalSize);
  sinfo.vector = _globalVec;
  sinfo.dm     = _dm;
  
  logger.stagePop();
} // createScatter

// ----------------------------------------------------------------------
// Create PETSc vector scatter for field. This is used to transfer
// information from the "global" PETSc vector view to the "local"
// Sieve section view. The PETSc vector does not contain constrained
// DOF. Use createScatterWithBC() to include the constrained DOF in
// the PETSc vector.
template<typename mesh_type>
template<typename scatter_mesh_type>
void
pylith::topology::Field<mesh_type>::createScatter(const scatter_mesh_type& mesh,
						  const typename ALE::Obj<typename SieveMesh::numbering_type> numbering,
						  const char* context)
{ // createScatter
  assert(!numbering.isNull());
  assert(context);
  PetscErrorCode err = 0;

  const bool createScatterOk = true;
  ScatterInfo& sinfo = _getScatter(context, createScatterOk);
  
  // Only create if scatter and scatterVec do not alreay exist.
  if (sinfo.scatter) {
    assert(sinfo.scatterVec);
    assert(sinfo.vector);
    return;
  } // if

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("GlobalOrder");
  PetscInt localSize, globalSize;

  err = PetscObjectReference((PetscObject) _dm);CHECK_PETSC_ERROR(err);
  err = PetscObjectReference((PetscObject) _globalVec);CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) _globalVec, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  err = VecGetSize(_localVec,  &localSize);CHECK_PETSC_ERROR(err);
  err = VecGetSize(_globalVec, &globalSize);CHECK_PETSC_ERROR(err);
  //assert(order->getLocalSize()  == localSize);
  //assert(order->getGlobalSize() == globalSize);
  sinfo.vector = _globalVec;
  sinfo.dm     = _dm;

#if 0
  std::cout << "CONTEXT: " << context 
	    << ", orderLabel: " << orderLabel
	    << ", section size w/BC: " << _section->sizeWithBC()
	    << ", section size: " << _section->size()
	    << ", global numbering size: " << numbering->getGlobalSize()
	    << ", global size: " << order->getGlobalSize()
	    << std::endl;
#endif
  
  logger.stagePop();
} // createScatter

// ----------------------------------------------------------------------
// Create PETSc vector scatter for field. This is used to transfer
// information from the "global" PETSc vector view to the "local"
// Sieve section view. The PETSc vector does not contain constrained
// DOF. Use createScatterWithBC() to include the constrained DOF in
// the PETSc vector.
template<typename mesh_type>
template<typename scatter_mesh_type>
void
pylith::topology::Field<mesh_type>::createScatterWithBC(const scatter_mesh_type& mesh,
							const char* context)
{ // createScatterWithBC
  assert(context);
  PetscErrorCode err = 0;

  const bool createScatterOk = true;
  ScatterInfo& sinfo = _getScatter(context, createScatterOk);
  if (sinfo.scatter) {
    assert(sinfo.scatterVec);
    assert(sinfo.vector);
    return;
  } // if

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("GlobalOrder");

  PetscSection section, newSection, gsection;
  PetscSF      sf;

  err = DMPlexClone(_dm, &sinfo.dm);CHECK_PETSC_ERROR(err);
  err = DMGetDefaultSection(_dm, &section);CHECK_PETSC_ERROR(err);
  err = PetscSectionClone(section, &newSection);CHECK_PETSC_ERROR(err);
  err = DMSetDefaultSection(sinfo.dm, newSection);CHECK_PETSC_ERROR(err);
  err = DMGetPointSF(sinfo.dm, &sf);CHECK_PETSC_ERROR(err);
  err = PetscSectionCreateGlobalSection(section, sf, PETSC_TRUE, &gsection);CHECK_PETSC_ERROR(err);
  err = DMSetDefaultGlobalSection(sinfo.dm, gsection);CHECK_PETSC_ERROR(err);
  err = DMCreateGlobalVector(sinfo.dm, &sinfo.vector);CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) sinfo.vector, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  PetscInt localSize, globalSize;

  err = PetscSectionGetStorageSize(section, &localSize);CHECK_PETSC_ERROR(err);
  err = VecGetSize(sinfo.vector, &globalSize);CHECK_PETSC_ERROR(err);
  //assert(order->getLocalSize()  == localSize);
  //assert(order->getGlobalSize() == globalSize);

  logger.stagePop();
} // createScatterWithBC

// ----------------------------------------------------------------------
// Create PETSc vector scatter for field. This is used to transfer
// information from the "global" PETSc vector view to the "local"
// Sieve section view. The PETSc vector includes constrained DOF. Use
// createScatter() if constrained DOF should be omitted from the PETSc
// vector.
template<typename mesh_type>
template<typename scatter_mesh_type>
void
pylith::topology::Field<mesh_type>::createScatterWithBC(const scatter_mesh_type& mesh,
							const std::string& labelName,
							PetscInt labelValue,
							const char* context)
{ // createScatterWithBC
  assert(context);
  PetscErrorCode err = 0;

  const bool createScatterOk = true;
  ScatterInfo& sinfo = _getScatter(context, createScatterOk);
  
  // Only create if scatter and scatterVec do not alreay exist.
  if (sinfo.scatter) {
    assert(sinfo.scatterVec);
    assert(sinfo.vector);
    return;
  } // if

  PetscDM dm = mesh.dmMesh();assert(dm);
  PetscSection section, newSection, gsection, subSection = NULL;
  PetscSF sf;
  PetscDMLabel subpointMap, subpointMapF;
  PetscInt dim, dimF, pStart, pEnd, qStart, qEnd, cEnd, cMax, vEnd, vMax;
  err = DMPlexGetHeightStratum(_dm, 0, NULL, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetDepthStratum(_dm, 0, NULL, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetHybridBounds(_dm, &cMax, NULL, NULL, &vMax);CHECK_PETSC_ERROR(err);
  PetscInt excludeRanges[4] = {cMax, cEnd, vMax, vEnd};
  PetscInt numExcludes = (cMax >= 0 ? 1 : 0) + (vMax >= 0 ? 1 : 0);

  err = DMGetDefaultSection(_dm, &section);CHECK_PETSC_ERROR(err);
  err = DMPlexGetDimension(dm,  &dim);CHECK_PETSC_ERROR(err);
  err = DMPlexGetDimension(_dm, &dimF);CHECK_PETSC_ERROR(err);
  err = DMPlexGetChart(dm,  &pStart, &pEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetChart(_dm, &qStart, &qEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetSubpointMap(dm,  &subpointMap);CHECK_PETSC_ERROR(err);
  err = DMPlexGetSubpointMap(_dm, &subpointMapF);CHECK_PETSC_ERROR(err);
  if (((dim != dimF) || ((pEnd-pStart) < (qEnd-qStart))) && subpointMap && !subpointMapF) {
    const PetscInt *ind;
    PetscIS subpointIS;
    PetscInt n, q;

    err = PetscPrintf(PETSC_COMM_SELF, "Making translation PetscSection\n");CHECK_PETSC_ERROR(err);
    err = PetscSectionGetChart(section, &qStart, &qEnd);CHECK_PETSC_ERROR(err);
    err = DMPlexCreateSubpointIS(dm, &subpointIS);CHECK_PETSC_ERROR(err);
    err = ISGetLocalSize(subpointIS, &n);CHECK_PETSC_ERROR(err);
    err = ISGetIndices(subpointIS, &ind);CHECK_PETSC_ERROR(err);
    err = PetscSectionCreate(mesh.comm(), &subSection);CHECK_PETSC_ERROR(err);
    err = PetscSectionSetChart(subSection, pStart, pEnd);CHECK_PETSC_ERROR(err);
    for(q = qStart; q < qEnd; ++q) {
      PetscInt dof, off, p;

      err = PetscSectionGetDof(section, q, &dof);CHECK_PETSC_ERROR(err);
      if (dof) {
        err = PetscFindInt(q, n, ind, &p);CHECK_PETSC_ERROR(err);
        if ((p >= pStart) && (p < pEnd)) {
          err = PetscSectionSetDof(subSection, p, dof);CHECK_PETSC_ERROR(err);
          err = PetscSectionGetOffset(section, q, &off);CHECK_PETSC_ERROR(err);
          err = PetscSectionSetOffset(subSection, p, off);CHECK_PETSC_ERROR(err);
        } // if
      } // if
    } // for
    err = ISRestoreIndices(subpointIS, &ind);CHECK_PETSC_ERROR(err);
    err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
    /* No need to setup section */
    err = PetscSectionView(subSection, PETSC_VIEWER_STDOUT_WORLD);CHECK_PETSC_ERROR(err);
    section = subSection;
    /* There are no excludes for surface meshes */
    numExcludes = 0;
  } // if

  err = DMPlexClone(_dm, &sinfo.dm);CHECK_PETSC_ERROR(err);
  err = PetscSectionClone(section, &newSection);CHECK_PETSC_ERROR(err);
  err = DMSetDefaultSection(sinfo.dm, newSection);CHECK_PETSC_ERROR(err);
  err = DMGetPointSF(sinfo.dm, &sf);CHECK_PETSC_ERROR(err);
  if (labelName.empty()) {
    err = PetscSectionCreateGlobalSectionCensored(section, sf, PETSC_TRUE, numExcludes, excludeRanges, &gsection);CHECK_PETSC_ERROR(err);
  } else {
    DMLabel label;

    err = DMPlexGetLabel(sinfo.dm, labelName.c_str(), &label);CHECK_PETSC_ERROR(err);
    err = PetscSectionCreateGlobalSectionLabel(section, sf, PETSC_TRUE, label, labelValue, &gsection);CHECK_PETSC_ERROR(err);
  } // if/else
  if (((dim != dimF) || ((pEnd-pStart) < (qEnd-qStart))) && subpointMap && !subpointMapF) {
    err = PetscSectionView(gsection, PETSC_VIEWER_STDOUT_WORLD);CHECK_PETSC_ERROR(err);
  } // if
  err = DMSetDefaultGlobalSection(sinfo.dm, gsection);CHECK_PETSC_ERROR(err);
  err = DMCreateGlobalVector(sinfo.dm, &sinfo.vector);CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) sinfo.vector, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  PetscInt localSize, globalSize;

  err = PetscSectionGetStorageSize(section, &localSize);CHECK_PETSC_ERROR(err);
  err = VecGetSize(sinfo.vector, &globalSize);CHECK_PETSC_ERROR(err);
  /* assert(order->getLocalSize()  == localSize); This does not work because the local vector includes the lagrange cell variables */
  /* assert(order->getGlobalSize() == globalSize); */
  if (subSection) {err = PetscSectionDestroy(&subSection);CHECK_PETSC_ERROR(err);}
#if 0
  std::cout << "["<<sieveMesh->commRank()<<"] CONTEXT: " << context 
	    << ", orderLabel: " << orderLabel
	    << ", section size w/BC: " << _section->sizeWithBC()
	    << ", section size: " << _section->size()
	    << ", section storage size: " << _section->getStorageSize()
	    << ", global numbering size: " << numbering->getGlobalSize()
	    << ", global order size: " << order->getGlobalSize()
	    << ", local numbering size: " << numbering->getLocalSize()
	    << ", local order size: " << order->getLocalSize()
	    << ", scatter from size: " << sinfo.scatter->from_n
	    << ", scatter: " << sinfo.scatter
	    << std::endl;
#endif
} // createScatterWithBC

// ----------------------------------------------------------------------
// Get PETSc vector associated with field.
template<typename mesh_type>
PetscVec
pylith::topology::Field<mesh_type>::vector(const char* context)
{ // vector
  std::ostringstream msg;

  ScatterInfo& sinfo = _getScatter(context);
  return sinfo.vector;
} // vector

// ----------------------------------------------------------------------
// Get PETSc vector associated with field.
template<typename mesh_type>
const PetscVec
pylith::topology::Field<mesh_type>::vector(const char* context) const
{ // vector
  std::ostringstream msg;

  const ScatterInfo& sinfo = _getScatter(context);
  return sinfo.vector;
} // vector

// ----------------------------------------------------------------------
// Scatter section information across processors to update the
//  PETSc vector view of the field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::scatterSectionToVector(const char* context) const
{ // scatterSectionToVector
  assert(context);

  const ScatterInfo& sinfo = _getScatter(context);
  scatterSectionToVector(sinfo.vector, context);
} // scatterSectionToVector

// ----------------------------------------------------------------------
// Scatter section information across processors to update the
//  PETSc vector view of the field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::scatterSectionToVector(const PetscVec vector,
							   const char* context) const
{ // scatterSectionToVector
  assert(vector);
  assert(context);
  const ScatterInfo& sinfo = _getScatter(context);
  PetscErrorCode err   = 0;
#if 0 // OBSOLETE??
  if (!_section.isNull()) {
    err = VecScatterBegin(sinfo.scatter, sinfo.scatterVec, vector,
                          INSERT_VALUES, SCATTER_FORWARD);CHECK_PETSC_ERROR(err);
    err = VecScatterEnd(sinfo.scatter, sinfo.scatterVec, vector,
                        INSERT_VALUES, SCATTER_FORWARD);CHECK_PETSC_ERROR(err);
  }
#endif
  if (sinfo.dm) {
    err = DMLocalToGlobalBegin(sinfo.dm, _localVec, INSERT_VALUES, vector);CHECK_PETSC_ERROR(err);
    err = DMLocalToGlobalEnd(sinfo.dm, _localVec, INSERT_VALUES, vector);CHECK_PETSC_ERROR(err);
  } // if
} // scatterSectionToVector

// ----------------------------------------------------------------------
// Scatter PETSc vector information across processors to update the
// section view of the field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::scatterVectorToSection(const char* context) const
{ // scatterVectorToSection
  assert(context);

  const ScatterInfo& sinfo = _getScatter(context);
  scatterVectorToSection(sinfo.vector, context);
} // scatterVectorToSection

// ----------------------------------------------------------------------
// Scatter PETSc vector information across processors to update the
// section view of the field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::scatterVectorToSection(const PetscVec vector,
									 const char* context) const
{ // scatterVectorToSection
  assert(vector);
  assert(context);
  const ScatterInfo& sinfo = _getScatter(context);
  PetscErrorCode err = 0;

#if 0 // OBSOLETE??
  if (!_section.isNull()) {
    err = VecScatterBegin(sinfo.scatter, vector, sinfo.scatterVec,
                          INSERT_VALUES, SCATTER_REVERSE); CHECK_PETSC_ERROR(err);
    err = VecScatterEnd(sinfo.scatter, vector, sinfo.scatterVec,
                        INSERT_VALUES, SCATTER_REVERSE); CHECK_PETSC_ERROR(err);
  }
#endif
  if (sinfo.dm) {
    err = DMGlobalToLocalBegin(sinfo.dm, vector, INSERT_VALUES, _localVec);CHECK_PETSC_ERROR(err);
    err = DMGlobalToLocalEnd(sinfo.dm, vector, INSERT_VALUES, _localVec);CHECK_PETSC_ERROR(err);
  } // if
} // scatterVectorToSection

// ----------------------------------------------------------------------
// Get fiber dimension associated with section (only works if fiber
// dimension is uniform).
template<typename mesh_type>
int
pylith::topology::Field<mesh_type>::_getFiberDim(void)
{ // _getFiberDim
  assert(_dm);

  PetscSection s = NULL;
  PetscInt pStart, pEnd;
  int fiberDimLocal, fiberDim = 0;
  PetscErrorCode err;

  err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetChart(s, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
  if (pEnd > pStart) {err = PetscSectionGetDof(s, pStart, &fiberDimLocal);CHECK_PETSC_ERROR(err);}
  MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPI_INT, MPI_MAX, _mesh.comm());

  return fiberDim;
} // _getFiberDim

// ----------------------------------------------------------------------
// Get scatter for given context.
template<typename mesh_type>
typename pylith::topology::Field<mesh_type>::ScatterInfo&
pylith::topology::Field<mesh_type>::_getScatter(const char* context,
						const bool createOk)
{ // _getScatter
  assert(context);

  bool isNewScatter = _scatters.find(context) == _scatters.end();

  // Synchronize creation of scatter (empty sections may have
  // leftover, reusable scatters that need to be cleared out).
  int numNewScatterLocal = (isNewScatter) ? 1 : 0;
  int numNewScatter = 0;
  MPI_Allreduce(&numNewScatterLocal, &numNewScatter, 1, MPI_INT, MPI_MAX, _mesh.comm());
  if (numNewScatter && !isNewScatter) {
    // remove old scatter
    ScatterInfo& sinfo = _scatters[context];
    PetscErrorCode err = 0;
    if (sinfo.vector) {
      err = VecDestroy(&sinfo.vector);CHECK_PETSC_ERROR(err);
    } // if
    if (sinfo.scatter) {
      err = VecScatterDestroy(&sinfo.scatter);CHECK_PETSC_ERROR(err);
    } // if

    if (sinfo.scatterVec) {
      err = VecDestroy(&sinfo.scatterVec);CHECK_PETSC_ERROR(err);
    } // if

    _scatters.erase(context);
    isNewScatter = true;
  } // if

  if (isNewScatter && !createOk) {
    std::ostringstream msg;
    msg << "Scatter for context '" << context << "' does not exist.";
    throw std::runtime_error(msg.str());
  } // if
  
  ScatterInfo& sinfo = _scatters[context];
  if (isNewScatter) {
    sinfo.dm = 0;
    sinfo.vector = 0;
    sinfo.scatter = 0;
    sinfo.scatterVec = 0;
  } // if
  assert(_scatters.find(context) != _scatters.end());

  return sinfo;
} // _getScatter

// ----------------------------------------------------------------------
// Get scatter for given context.
template<typename mesh_type>
const typename pylith::topology::Field<mesh_type>::ScatterInfo&
pylith::topology::Field<mesh_type>::_getScatter(const char* context) const
{ // _getScatter
  assert(context);

  const typename scatter_map_type::const_iterator s_iter = 
    _scatters.find(context);
  if (s_iter == _scatters.end()) {
    std::ostringstream msg;
    msg << "Scatter for context '" << context << "' does not exist.";
    throw std::runtime_error(msg.str());
  } // if
  
  return s_iter->second;
} // _getScatter

// ----------------------------------------------------------------------
// Experimental
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::addField(const char *name,
					     int numComponents)
{
  // Keep track of name/components until setup
  _tmpFields[name] = numComponents;
  _metadata[name]  = _metadata["default"];
}

// ----------------------------------------------------------------------
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::setupFields()
{
  assert(_dm);
  // Keep track of name/components until setup
  PetscSection section;
  PetscInt f = 0;
  PetscErrorCode err = DMGetDefaultSection(_dm, &section);CHECK_PETSC_ERROR(err);assert(section);
  err = PetscSectionSetNumFields(section, _tmpFields.size());CHECK_PETSC_ERROR(err);
  for(std::map<std::string, int>::const_iterator f_iter = _tmpFields.begin(); f_iter != _tmpFields.end(); ++f_iter, ++f) {
    err = PetscSectionSetFieldName(section, f, f_iter->first.c_str());CHECK_PETSC_ERROR(err);
    err = PetscSectionSetFieldComponents(section, f, f_iter->second);CHECK_PETSC_ERROR(err);
  } // for
  _tmpFields.clear();
#if 0 // :MATT: What is going on here? Is this obsolete?
  // Right now, we assume that the section covers the entire chart
  PetscInt pStart, pEnd;

  err = DMPlexGetChart(_dm, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
  err = PetscSectionSetChart(section, pStart, pEnd);CHECK_PETSC_ERROR(err);
#endif
}

// ----------------------------------------------------------------------
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::updateDof(const char *name,
					      const DomainEnum domain,
					      int fiberDim)
{ // updateDof
  PetscInt pStart, pEnd, f = 0;
  PetscErrorCode err;

  assert(_dm);
  switch(domain) {
  case VERTICES_FIELD:
    err = DMPlexGetDepthStratum(_dm, 0, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    break;
  case CELLS_FIELD:
    err = DMPlexGetHeightStratum(_dm, 0, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    break;
  case FACES_FIELD:
    err = DMPlexGetHeightStratum(_dm, 1, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    break;
  case POINTS_FIELD:
    err = DMPlexGetChart(_dm, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    break;
  default:
    std::cerr << "Unknown value for DomainEnum: " << domain << std::endl;
    throw std::logic_error("Bad domain enum in Field.");
  }
  PetscSection section = NULL;
  err = DMGetDefaultSection(_dm, &section);CHECK_PETSC_ERROR(err);assert(section);
  for(map_type::const_iterator f_iter = _metadata.begin(); f_iter != _metadata.end(); ++f_iter) {
    if (f_iter->first == name) break;
    if (f_iter->first == "default") continue;
    ++f;
  } // for
  assert(f < _metadata.size());
  for(PetscInt p = pStart; p < pEnd; ++p) {
    //err = PetscSectionAddDof(section, p, fiberDim);CHECK_PETSC_ERROR(err);
    err = PetscSectionSetFieldDof(section, p, f, fiberDim);CHECK_PETSC_ERROR(err);
  } // for
} // updateDof

// End of file 
