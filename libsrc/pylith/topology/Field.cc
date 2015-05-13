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

#include "Field.hh" // implementation of class methods

#include "Mesh.hh" // USES Mesh

#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR
#include <iostream> // USES std::cout

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::Field::Field(const Mesh& mesh) :
  _mesh(mesh),
  _dm(NULL),
  _globalVec(NULL),
  _localVec(NULL)
{ // constructor
  PYLITH_METHOD_BEGIN;

  _metadata.label = "unknown";
  _metadata.vectorFieldType = OTHER;
  _metadata.scale = 1.0;
  _metadata.dimsOkay = false;
  if (mesh.dmMesh()) {
    PetscDM dm = mesh.dmMesh();assert(dm);
    PetscVec coordVec = NULL;
    PetscSection s = NULL;
    PetscErrorCode err;

    err = DMDestroy(&_dm);PYLITH_CHECK_ERROR(err);
    err = DMClone(dm, &_dm);PYLITH_CHECK_ERROR(err);
    err = DMGetCoordinatesLocal(dm, &coordVec);PYLITH_CHECK_ERROR(err);
    if (coordVec) {
      PetscDM coordDM=NULL, newCoordDM=NULL;
      PetscSection coordSection=NULL, newCoordSection=NULL;

      err = DMGetCoordinateDM(dm, &coordDM);PYLITH_CHECK_ERROR(err);
      err = DMGetCoordinateDM(_dm, &newCoordDM);PYLITH_CHECK_ERROR(err);
      err = DMGetDefaultSection(coordDM, &coordSection);PYLITH_CHECK_ERROR(err);
      err = PetscSectionClone(coordSection, &newCoordSection);PYLITH_CHECK_ERROR(err);
      err = DMSetDefaultSection(newCoordDM, newCoordSection);PYLITH_CHECK_ERROR(err);
      err = PetscSectionDestroy(&newCoordSection);PYLITH_CHECK_ERROR(err);
      err = DMSetCoordinatesLocal(_dm, coordVec);PYLITH_CHECK_ERROR(err);
    } // if
    err = PetscSectionCreate(mesh.comm(), &s);PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultSection(_dm, s);PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&s);PYLITH_CHECK_ERROR(err);
  } // if

  PYLITH_METHOD_END;
} // constructor

// ----------------------------------------------------------------------
// Constructor with mesh, DM, and metadata
pylith::topology::Field::Field(const Mesh& mesh,
			       PetscDM dm,
			       const Metadata& metadata) :
  _mesh(mesh),
  _dm(dm),
  _globalVec(NULL),
  _localVec(NULL)
{ // constructor
  PYLITH_METHOD_BEGIN;

  assert(dm);
  PetscErrorCode err;

  _metadata = metadata;
  err = DMCreateLocalVector(_dm, &_localVec);PYLITH_CHECK_ERROR(err);
  err = DMCreateGlobalVector(_dm, &_globalVec);PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetName((PetscObject) _globalVec, _metadata.label.c_str());PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetName((PetscObject) _localVec,  _metadata.label.c_str());PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // constructor

// ----------------------------------------------------------------------
// Constructor with mesh, DM, local data, and metadata
pylith::topology::Field::Field(const Mesh& mesh,
			       PetscDM dm,
			       PetscVec localVec,
			       const Metadata& metadata) :
  _mesh(mesh),
  _dm(dm),
  _globalVec(NULL),
  _localVec(NULL)
{ // constructor
  PYLITH_METHOD_BEGIN;

  assert(dm);
  assert(localVec);

  PetscErrorCode err;

  _metadata = metadata;
  err = DMCreateLocalVector(_dm, &_localVec);PYLITH_CHECK_ERROR(err);
  err = DMCreateGlobalVector(_dm, &_globalVec);PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetName((PetscObject) _globalVec, _metadata.label.c_str());PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetName((PetscObject) _localVec,  _metadata.label.c_str());PYLITH_CHECK_ERROR(err);
  err = VecCopy(localVec, _localVec);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::Field::~Field(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::Field::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  clear();
  PetscErrorCode err = DMDestroy(&_dm);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set label for field.
void
pylith::topology::Field::label(const char* value)
{ // label
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err;

  _metadata.label = value;
  if (_localVec)  {
    err = PetscObjectSetName((PetscObject) _localVec, value);PYLITH_CHECK_ERROR(err);
  } // if
  if (_globalVec) {
    err = PetscObjectSetName((PetscObject) _globalVec, value);PYLITH_CHECK_ERROR(err);
  } // if

  const scatter_map_type::const_iterator scattersEnd = _scatters.end();
  for (scatter_map_type::const_iterator s_iter=_scatters.begin();
       s_iter != scattersEnd;
       ++s_iter) {
    if (s_iter->second.vector) {
      err = PetscObjectSetName((PetscObject)s_iter->second.vector, value);PYLITH_CHECK_ERROR(err);    
    } // if
  } // for

  PYLITH_METHOD_END;
} // label

// ----------------------------------------------------------------------
// Get spatial dimension of domain.
int
pylith::topology::Field::spaceDim(void) const
{ // spaceDim
  const spatialdata::geocoords::CoordSys* cs = _mesh.coordsys();
  return (cs) ? cs->spaceDim() : 0;
} // spaceDim

// ----------------------------------------------------------------------
// Has section been setup?
bool
pylith::topology::Field::hasSection(void) const
{ // hasSection
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err;
  PetscSection s = NULL;
  err = DMGetDefaultSection(_dm, &s);PYLITH_CHECK_ERROR(err);

  PetscInt pStart = 0, pEnd = 0;
  err = PetscSectionGetChart(s, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
  
  bool result = (pEnd < 0) ? false : true;
  
  PYLITH_METHOD_RETURN(result);
} // hasSection

// ----------------------------------------------------------------------
// Get the chart size.
int
pylith::topology::Field::chartSize(void) const
{ // chartSize
  PYLITH_METHOD_BEGIN;

  assert(_dm);
  PetscSection s = NULL;
  PetscInt pStart, pEnd;
  PetscErrorCode err;

  err = DMGetDefaultSection(_dm, &s);PYLITH_CHECK_ERROR(err);
  err = PetscSectionGetChart(s, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_RETURN(pEnd-pStart);
} // chartSize

// ----------------------------------------------------------------------
// Get the number of degrees of freedom.
int
pylith::topology::Field::sectionSize(void) const
{ // sectionSize
  PYLITH_METHOD_BEGIN;

  PetscInt size = 0;

  if (_dm) {
    PetscSection s = NULL;
    PetscErrorCode err;

    err = DMGetDefaultSection(_dm, &s);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetStorageSize(s, &size);PYLITH_CHECK_ERROR(err);
  } // if

  PYLITH_METHOD_RETURN(size);
} // sectionSize

// ----------------------------------------------------------------------
// Set chart for solution.
void
pylith::topology::Field::setupSolnChart(void)
{ // setupSolnChart
  PYLITH_METHOD_BEGIN;

  assert(_dm);

  // :TODO: Update this to use discretization information after removing FIAT.

  // :KLUDGE: Assume solution has DOF over vertices and hybrid edges.
  PetscErrorCode err;
  // Get range of vertices.
  PetscInt pStart = -1;
  PetscInt pEnd = -1;
  err = DMPlexGetDepthStratum(_dm, 0, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
  // Get last edge and 
  PetscInt eEnd = -1;
  PetscInt eMax = -1;
  err = DMPlexGetDepthStratum(_dm, 1, NULL, &eEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetHybridBounds(_dm, NULL, NULL, &eMax, NULL);PYLITH_CHECK_ERROR(err);
  // If have hybrid edges, extend chart to include hybrid edges; otherwise just use points.
  if (eEnd > eMax) {
    pEnd = eEnd;
  } // if

  PetscSection s = NULL;
  if (pStart < pEnd) {
    err = DMGetDefaultSection(_dm, &s);PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultGlobalSection(_dm, NULL);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetChart(s, pStart, pEnd);PYLITH_CHECK_ERROR(err);
  } else { // create empty chart
    err = DMGetDefaultSection(_dm, &s);PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultGlobalSection(_dm, NULL);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetChart(s, 0, 0);PYLITH_CHECK_ERROR(err);
  } // if/else

  PYLITH_METHOD_END;
} // setupSolnChart


// ----------------------------------------------------------------------
// Setup default DOF for solution.
void
  pylith::topology::Field::setupSolnDof(const int fiberDim,
					const char* subfieldName)
{ // setupSolnDof
  PYLITH_METHOD_BEGIN;

  assert(_dm);

  // :TODO: Update this to use discretization information after removing FIAT.

  // :KLUDGE: Assume solution has DOF over vertices.
  const int subfieldIndex = _subfields[subfieldName].index;

  PetscErrorCode err;
  // Get range of vertices.
  PetscInt pStart = -1;
  PetscInt pEnd = -1;
  err = DMPlexGetDepthStratum(_dm, 0, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

  PetscSection s = NULL;
  err = DMGetDefaultSection(_dm, &s);PYLITH_CHECK_ERROR(err);
  for (PetscInt p = pStart; p < pEnd; ++p) {
    err = PetscSectionSetDof(s, p, fiberDim);PYLITH_CHECK_ERROR(err);

    // Set DOF in subfield
    err = PetscSectionSetFieldDof(s, p, subfieldIndex, fiberDim);PYLITH_CHECK_ERROR(err);
  } // for

  PYLITH_METHOD_END;
} // setupSolnDof


// ----------------------------------------------------------------------
// Create PETSc section and set chart and fiber dimesion for a list of
// points.
void
pylith::topology::Field::newSection(const int_array& points,
				    const int fiberDim)
{ // newSection
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err;

  // Clear memory
  clear();
  assert(_dm);
  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _metadata.label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if
  
  const PetscInt npts = points.size();
  if (npts > 0) {
    PetscSection s = NULL;
    PetscInt pointMin = points[0], pointMax = points[0];

    for (PetscInt i = 1; i < npts; ++i) {
      pointMin = std::min(pointMin, points[i]);
      pointMax = std::max(pointMax, points[i]);
    } // for
    err = DMGetDefaultSection(_dm, &s);PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultGlobalSection(_dm, NULL);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetChart(s, pointMin, pointMax+1);PYLITH_CHECK_ERROR(err);
    for (PetscInt i = 0; i < npts; ++i) {
      err = PetscSectionSetDof(s, points[i], fiberDim);PYLITH_CHECK_ERROR(err);
    } // for
  } else { // create empty chart
    PetscSection s = NULL;

    err = DMGetDefaultSection(_dm, &s);PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultGlobalSection(_dm, NULL);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetChart(s, 0, 0);PYLITH_CHECK_ERROR(err);
  } // if/else

  PYLITH_METHOD_END;
} // newSection

// ----------------------------------------------------------------------
// Create PETSc section and set chart and fiber dimesion for a list of
// points.
void
pylith::topology::Field::newSection(const PetscInt *points, 
				    const PetscInt num,
				    const int fiberDim)
{ // newSection
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err;

  // Clear memory
  clear();
  assert(_dm);
  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _metadata.label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if
  
  if (num > 0) {
    PetscSection s = NULL;
    PetscInt pointMin = points[0], pointMax = points[0];

    for (PetscInt i = 1; i < num; ++i) {
      pointMin = std::min(pointMin, points[i]);
      pointMax = std::max(pointMax, points[i]);
    } // for
    err = DMGetDefaultSection(_dm, &s);PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultGlobalSection(_dm, NULL);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetChart(s, pointMin, pointMax+1);PYLITH_CHECK_ERROR(err);
    for (PetscInt i = 0; i < num; ++i) {
      err = PetscSectionSetDof(s, points[i], fiberDim);PYLITH_CHECK_ERROR(err);
    } // for
  } else { // create empty chart
    PetscSection s = NULL;

    err = DMGetDefaultSection(_dm, &s);PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultGlobalSection(_dm, NULL);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetChart(s, 0, 0);PYLITH_CHECK_ERROR(err);
  } // if/else

  PYLITH_METHOD_END;
} // newSection

// ----------------------------------------------------------------------
// Create PETSc section and set chart and fiber dimesion.
void
pylith::topology::Field::newSection(const DomainEnum domain,
				    const int fiberDim,
				    const int stratum)
{ // newSection
  PYLITH_METHOD_BEGIN;

  // Changing this because cells/vertices are numbered differently in the new scheme
  assert(_dm);
  PetscInt pStart, pEnd;
  PetscErrorCode err;

  switch(domain) {
  case VERTICES_FIELD:
    err = DMPlexGetDepthStratum(_dm, stratum, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    break;
  case CELLS_FIELD:
    err = DMPlexGetHeightStratum(_dm, stratum, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    break;
  case FACES_FIELD:
    err = DMPlexGetHeightStratum(_dm, stratum+1, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    break;
  case POINTS_FIELD:
    err = DMPlexGetChart(_dm, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    break;
  default:
    std::ostringstream msg;
    msg << "Unknown value for DomainEnum: " << domain << "  in Field" << std::endl;
    throw std::logic_error(msg.str());
  }
  newSection(pStart, pEnd, fiberDim);

  PYLITH_METHOD_END;
} // newSection

// ----------------------------------------------------------------------
// Create PETSc section and set chart and fiber dimesion.
void
pylith::topology::Field::newSection(const PetscInt pStart, 
				    const PetscInt pEnd,
				    const int fiberDim)
{ // newSection
  PYLITH_METHOD_BEGIN;

  // Changing this because cells/vertices are numbered differently in the new scheme
  assert(_dm);
  PetscSection s = NULL;
  PetscErrorCode err;

  err = DMGetDefaultSection(_dm, &s);PYLITH_CHECK_ERROR(err);
  err = DMSetDefaultGlobalSection(_dm, NULL);PYLITH_CHECK_ERROR(err);
  err = PetscSectionSetChart(s, pStart, pEnd);PYLITH_CHECK_ERROR(err);

  for(PetscInt p = pStart; p < pEnd; ++p) {
    err = PetscSectionSetDof(s, p, fiberDim);PYLITH_CHECK_ERROR(err);
  } // for

  PYLITH_METHOD_END;
} // newSection

// ----------------------------------------------------------------------
// Create section given chart.
void
pylith::topology::Field::newSection(const Field& src,
				    const int fiberDim)
{ // newSection
  PYLITH_METHOD_BEGIN;

  // Clear memory
  clear();
  assert(_dm);assert(src._dm);

  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _metadata.label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  PetscSection srcs=NULL, s=NULL;
  PetscInt pStart, pEnd;
  PetscErrorCode err;

  err = DMSetDefaultGlobalSection(_dm, NULL);PYLITH_CHECK_ERROR(err);

  err = DMGetDefaultSection(_dm, &s);PYLITH_CHECK_ERROR(err);
  err = DMGetDefaultSection(src._dm, &srcs);PYLITH_CHECK_ERROR(err);
  err = PetscSectionGetChart(srcs, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
  err = PetscSectionSetChart(s, pStart, pEnd);PYLITH_CHECK_ERROR(err);
  for(PetscInt p = pStart; p < pEnd; ++p) {
    err = PetscSectionSetDof(s, p, fiberDim);PYLITH_CHECK_ERROR(err);
  } // for

  PYLITH_METHOD_END;
} // newSection

// ----------------------------------------------------------------------
// Create section with same layout as another section.
void
pylith::topology::Field::cloneSection(const Field& src)
{ // cloneSection
  PYLITH_METHOD_BEGIN;

  std::string origLabel = _metadata.label;

  // Clear memory
  clear();

  _metadata = src._metadata;
  label(origLabel.c_str());

  PetscSection section = src.localSection();
  PetscSection newSection = NULL;
  PetscErrorCode err;

  assert(_dm);
  err = PetscSectionClone(section, &newSection);PYLITH_CHECK_ERROR(err);
  err = DMSetDefaultSection(_dm, newSection);PYLITH_CHECK_ERROR(err);
  err = PetscSectionDestroy(&newSection);PYLITH_CHECK_ERROR(err);

  assert(!_globalVec);
  err = DMCreateGlobalVector(_dm, &_globalVec);PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetName((PetscObject) _globalVec, _metadata.label.c_str());PYLITH_CHECK_ERROR(err);

  assert(!_localVec);
  err = DMCreateLocalVector(_dm, &_localVec);PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetName((PetscObject) _localVec,  _metadata.label.c_str());PYLITH_CHECK_ERROR(err);
    
  // Reuse scatters in clone
  _scatters.clear();
  const scatter_map_type::const_iterator scattersEnd = src._scatters.end();
  for (scatter_map_type::const_iterator s_iter=src._scatters.begin(); s_iter != scattersEnd; ++s_iter) {
    ScatterInfo& sinfo = _scatters[s_iter->first];
    sinfo.dm = 0;
    sinfo.vector = 0;

    // Copy DM
    sinfo.dm = s_iter->second.dm;
    err = PetscObjectReference((PetscObject) sinfo.dm);PYLITH_CHECK_ERROR(err);

    // Create vector using sizes from source section
    PetscInt vecGlobalSize = 0, vecGlobalSize2 = 0;
    err = VecGetSize(s_iter->second.vector, &vecGlobalSize);PYLITH_CHECK_ERROR(err);
    err = VecGetSize(_globalVec, &vecGlobalSize2);PYLITH_CHECK_ERROR(err);      
    if (vecGlobalSize != vecGlobalSize2) {
      err = DMCreateGlobalVector(sinfo.dm, &sinfo.vector);PYLITH_CHECK_ERROR(err);
    } else {
      sinfo.vector = _globalVec;
      err = PetscObjectReference((PetscObject) sinfo.vector);PYLITH_CHECK_ERROR(err);
    } // if/else
    err = PetscObjectSetName((PetscObject)sinfo.vector, _metadata.label.c_str());PYLITH_CHECK_ERROR(err);
  } // for

  // Reuse subfields in clone
  _subfields.clear();
  const subfields_type::const_iterator subfieldsEnd = src._subfields.end();
  for (subfields_type::const_iterator s_iter=src._subfields.begin(); s_iter != subfieldsEnd; ++s_iter) {
    SubfieldInfo& sinfo = _subfields[s_iter->first];
    sinfo.metadata = s_iter->second.metadata;
    sinfo.index = s_iter->second.index;
    sinfo.dm = s_iter->second.dm;
    if (sinfo.dm) {
      err = PetscObjectReference((PetscObject) sinfo.dm);PYLITH_CHECK_ERROR(err);
    } // if
  } // for

  PYLITH_METHOD_END;
} // cloneSection

// ----------------------------------------------------------------------
// Clear variables associated with section.
void
pylith::topology::Field::clear(void)
{ // clear
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err = 0;
  
  const scatter_map_type::const_iterator scattersEnd = _scatters.end();
  for (scatter_map_type::iterator s_iter=_scatters.begin(); s_iter != scattersEnd; ++s_iter) {
    err = DMDestroy(&s_iter->second.dm);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&s_iter->second.vector);PYLITH_CHECK_ERROR(err);
  } // for
  _scatters.clear();

  err = VecDestroy(&_globalVec);PYLITH_CHECK_ERROR(err);
  err = VecDestroy(&_localVec);PYLITH_CHECK_ERROR(err);

  const subfields_type::const_iterator subfieldsEnd = _subfields.end();
  for (subfields_type::iterator s_iter=_subfields.begin(); s_iter != subfieldsEnd; ++s_iter) {
    err = DMDestroy(&s_iter->second.dm);PYLITH_CHECK_ERROR(err);
  } // for

  _metadata.scale = 1.0;
  _metadata.vectorFieldType = OTHER;
  _metadata.dimsOkay = false;

  PYLITH_METHOD_END;
} // clear

// ----------------------------------------------------------------------
// Allocate PETSc section.
void
pylith::topology::Field::allocate(void)
{ // allocate
  PYLITH_METHOD_BEGIN;

  PetscSection s = NULL;
  PetscErrorCode err;

  if (_dm) {
    err = DMGetDefaultSection(_dm, &s);PYLITH_CHECK_ERROR(err);
  } // if
  assert(s);
  err = PetscSectionSetUp(s);PYLITH_CHECK_ERROR(err);

  err = VecDestroy(&_globalVec);PYLITH_CHECK_ERROR(err);
  err = DMCreateGlobalVector(_dm, &_globalVec);PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetName((PetscObject) _globalVec, _metadata.label.c_str());PYLITH_CHECK_ERROR(err);

  err = VecDestroy(&_localVec);PYLITH_CHECK_ERROR(err);
  err = DMCreateLocalVector(_dm, &_localVec);PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetName((PetscObject) _localVec,  _metadata.label.c_str());PYLITH_CHECK_ERROR(err);

  // Create DM for subfields.
  int fields[1];
  for (subfields_type::iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
    SubfieldInfo& info = s_iter->second;
    fields[0] = info.index;
    err = DMDestroy(&info.dm);PYLITH_CHECK_ERROR(err);
    err = DMCreateSubDM(_dm, 1, fields, NULL, &info.dm);PYLITH_CHECK_ERROR(err);
  } // for
  

  PYLITH_METHOD_END;
} // allocate

// ----------------------------------------------------------------------
// Zero section values (excluding constrained DOF).
void
pylith::topology::Field::zero(void)
{ // zero
  PYLITH_METHOD_BEGIN;

  assert(_localVec);
  PetscSection section = NULL;
  PetscInt pStart, pEnd, maxDof = 0;
  PetscErrorCode err;

  err = DMGetDefaultSection(_dm, &section);PYLITH_CHECK_ERROR(err);    
  err = PetscSectionGetChart(section, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);    
  if (pEnd > pStart) {err = PetscSectionGetMaxDof(section, &maxDof);PYLITH_CHECK_ERROR(err);}
  scalar_array values(maxDof);
  values *= 0.0;

  for(PetscInt p = pStart; p < pEnd; ++p) {
    PetscInt dof;

    err = PetscSectionGetDof(section, p, &dof);PYLITH_CHECK_ERROR(err);
    if (dof > 0) {
      assert(dof <= maxDof);
      err = VecSetValuesSection(_localVec, section, p, &values[0], INSERT_VALUES);PYLITH_CHECK_ERROR(err);
    } // if
  } // for

  PYLITH_METHOD_END;
} // zero

// ----------------------------------------------------------------------
// Zero section values (including constrained DOF).
void
pylith::topology::Field::zeroAll(void)
{ // zeroAll
  PYLITH_METHOD_BEGIN;

  assert(_localVec);
  PetscErrorCode err = VecSet(_localVec, 0.0);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // zeroAll

// ----------------------------------------------------------------------
// Complete section by assembling across processors.
void
pylith::topology::Field::complete(void)
{ // complete
  PYLITH_METHOD_BEGIN;

  assert(_dm);
  // Not sure if DMLocalToLocal() would work
  PetscErrorCode err;

  err = VecSet(_globalVec, 0.0);PYLITH_CHECK_ERROR(err);
  err = DMLocalToGlobalBegin(_dm, _localVec, ADD_VALUES, _globalVec);PYLITH_CHECK_ERROR(err);
  err = DMLocalToGlobalEnd(_dm, _localVec, ADD_VALUES, _globalVec);PYLITH_CHECK_ERROR(err);
  err = DMGlobalToLocalBegin(_dm, _globalVec, INSERT_VALUES, _localVec);PYLITH_CHECK_ERROR(err);
  err = DMGlobalToLocalEnd(_dm, _globalVec, INSERT_VALUES, _localVec);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // complete

// ----------------------------------------------------------------------
// Copy field values and metadata.
void
pylith::topology::Field::copy(const Field& field)
{ // copy
  PYLITH_METHOD_BEGIN;

  // Check compatibility of sections
  const int srcSize = field.chartSize();
  const int dstSize = chartSize();
  if (field.spaceDim() != spaceDim() ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot copy values from section '" << _metadata.label 
	<< "' to section '" << _metadata.label
	<< "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    space dim: " << field.spaceDim() << "\n"
	<< "    vector field type: " << field._metadata.vectorFieldType << "\n"
	<< "    scale: " << field._metadata.scale << "\n"
	<< "    size: " << srcSize << "\n"
	<< "  Destination section:\n"
	<< "    space dim: " << spaceDim() << "\n"
	<< "    vector field type: " << _metadata.vectorFieldType << "\n"
	<< "    scale: " << _metadata.scale << "\n"
	<< "    size: " << dstSize;
    throw std::runtime_error(msg.str());
  } // if
  assert(_localVec && field._localVec);

  PetscErrorCode err = VecCopy(field._localVec, _localVec);PYLITH_CHECK_ERROR(err);

  // Update metadata
  label(field._metadata.label.c_str());
  _metadata.vectorFieldType = field._metadata.vectorFieldType;
  _metadata.scale = field._metadata.scale;

  PYLITH_METHOD_END;
} // copy

// ----------------------------------------------------------------------
// Add two fields, storing the result in one of the fields.
pylith::topology::Field&
pylith::topology::Field::operator+=(const Field& field)
{ // operator+=
  PYLITH_METHOD_BEGIN;

  // Check compatibility of sections
  const int srcSize = field.chartSize();
  const int dstSize = chartSize();
  if (field.spaceDim() != spaceDim() ||
      field._metadata.vectorFieldType != _metadata.vectorFieldType ||
      field._metadata.scale != _metadata.scale ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot add values from section '" << field._metadata.label 
	<< "' to section '" << _metadata.label
	<< "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    space dim: " << field.spaceDim() << "\n"
	<< "    vector field type: " << field._metadata.vectorFieldType << "\n"
    << "    scale: " << field._metadata.scale << "\n"
	<< "    size: " << srcSize << "\n"
	<< "  Destination section:\n"
	<< "    space dim: " << spaceDim() << "\n"
	<< "    vector field type: " << _metadata.vectorFieldType << "\n"
	<< "    scale: " << _metadata.scale << "\n"
	<< "    size: " << dstSize;
    throw std::runtime_error(msg.str());
  } // if
  assert(_localVec && field._localVec);
  PetscErrorCode err = VecAXPY(_localVec, 1.0, field._localVec);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_RETURN(*this);
} // operator+=

// ----------------------------------------------------------------------
// Dimensionalize field.
void
pylith::topology::Field::dimensionalize(void) const
{ // dimensionalize
  PYLITH_METHOD_BEGIN;

  if (!_metadata.dimsOkay) {
    std::ostringstream msg;
    msg << "Cannot dimensionalize field '" << _metadata.label
	<< "' because the flag has been set to keep field nondimensional.";
    throw std::runtime_error(msg.str());
  } // if

  spatialdata::units::Nondimensional normalizer;
  assert(_localVec);
  PetscSection section = NULL;
  PetscScalar *array = NULL;
  PetscInt pStart, pEnd;
  PetscErrorCode err;

  err = DMGetDefaultSection(_dm, &section);PYLITH_CHECK_ERROR(err);
  err = PetscSectionGetChart(section, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
  err = VecGetArray(_localVec, &array);PYLITH_CHECK_ERROR(err);
  for(PetscInt p = pStart; p < pEnd; ++p) {
    PetscInt dof, off;

    err = PetscSectionGetDof(section, p, &dof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(section, p, &off);PYLITH_CHECK_ERROR(err);
    if (dof) {
      normalizer.dimensionalize(&array[off], dof, _metadata.scale);
    }
  }
  err = VecRestoreArray(_localVec, &array);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // dimensionalize

// ----------------------------------------------------------------------
// Print field to standard out.
void
pylith::topology::Field::view(const char* label) const
{ // view
  PYLITH_METHOD_BEGIN;

  std::string vecFieldString;
  switch(_metadata.vectorFieldType)
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
    default:
      std::ostringstream msg;
      msg << "Unknown vector field value '" << _metadata.vectorFieldType << "'  in Field." << std::endl;
      throw std::logic_error(msg.str());
    } // switch

  std::cout << "Viewing field '" << _metadata.label << "' "<< label << ".\n"
	    << "  vector field type: " << vecFieldString << "\n"
	    << "  scale: " << _metadata.scale << "\n"
	    << "  dimensionalize flag: " << _metadata.dimsOkay << std::endl;
  if (_dm) {
    PetscSection   section = NULL;
    PetscMPIInt    numProcs, rank;
    PetscErrorCode err;

    err = DMGetDefaultSection(_dm, &section);PYLITH_CHECK_ERROR(err);
    err = DMView(_dm, PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
    err = PetscSectionView(section, PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
    err = MPI_Comm_size(PetscObjectComm((PetscObject) _dm), &numProcs);PYLITH_CHECK_ERROR(err);
    err = MPI_Comm_rank(PetscObjectComm((PetscObject) _dm), &rank);PYLITH_CHECK_ERROR(err);
    for (PetscInt p = 0; p < numProcs; ++p) {
      err = PetscPrintf(PetscObjectComm((PetscObject) _dm), "Proc %d local vector\n", p);PYLITH_CHECK_ERROR(err);
      if (p == rank) {err = VecView(_localVec, PETSC_VIEWER_STDOUT_SELF);PYLITH_CHECK_ERROR(err);}
      err = PetscBarrier((PetscObject) _dm);PYLITH_CHECK_ERROR(err);
    }
  }

  PYLITH_METHOD_END;
} // view

// ----------------------------------------------------------------------
// Create PETSc vector scatter for field. This is used to transfer
// information from the "global" PETSc vector view to the "local"
// PETSc section view.
void
pylith::topology::Field::createScatter(const Mesh& mesh,
				       const char* context)
{ // createScatter
  PYLITH_METHOD_BEGIN;

  assert(context);
  PetscErrorCode err = 0;

  const bool createScatterOk = true;
  ScatterInfo& sinfo = _getScatter(context, createScatterOk);
  if (sinfo.dm) {
    assert(sinfo.vector);
    PYLITH_METHOD_END;
  } // if

  err = DMDestroy(&sinfo.dm);PYLITH_CHECK_ERROR(err);
  sinfo.dm = _dm;
  err = PetscObjectReference((PetscObject) sinfo.dm);PYLITH_CHECK_ERROR(err);

  err = VecDestroy(&sinfo.vector);PYLITH_CHECK_ERROR(err);
  sinfo.vector = _globalVec;
  err = PetscObjectReference((PetscObject) sinfo.vector);PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetName((PetscObject) sinfo.vector, _metadata.label.c_str());PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // createScatter

// ----------------------------------------------------------------------
// Create PETSc vector scatter for field. This is used to transfer
// information from the "global" PETSc vector view to the "local"
// PETSc section view. The PETSc vector does not contain constrained
// DOF. Use createScatterWithBC() to include the constrained DOF in
// the PETSc vector.
void
pylith::topology::Field::createScatterWithBC(const Mesh& mesh,
					     const char* context)
{ // createScatterWithBC
  PYLITH_METHOD_BEGIN;

  assert(context);
  PetscErrorCode err = 0;

  const bool createScatterOk = true;
  ScatterInfo& sinfo = _getScatter(context, createScatterOk);
  if (sinfo.dm) {
    assert(sinfo.vector);
    PYLITH_METHOD_END;
  } // if

  PetscSection section = NULL, newSection = NULL, gsection = NULL;
  PetscSF sf = NULL;

  err = DMDestroy(&sinfo.dm);PYLITH_CHECK_ERROR(err);
  err = DMClone(_dm, &sinfo.dm);PYLITH_CHECK_ERROR(err);
  err = DMGetDefaultSection(_dm, &section);PYLITH_CHECK_ERROR(err);
  err = PetscSectionClone(section, &newSection);PYLITH_CHECK_ERROR(err);
  err = DMSetDefaultSection(sinfo.dm, newSection);PYLITH_CHECK_ERROR(err);
  err = PetscSectionDestroy(&newSection);PYLITH_CHECK_ERROR(err);
  err = DMGetPointSF(sinfo.dm, &sf);PYLITH_CHECK_ERROR(err);
  err = PetscSectionCreateGlobalSection(section, sf, PETSC_TRUE, &gsection);PYLITH_CHECK_ERROR(err);
  err = DMSetDefaultGlobalSection(sinfo.dm, gsection);PYLITH_CHECK_ERROR(err);
  err = PetscSectionDestroy(&gsection);PYLITH_CHECK_ERROR(err);
  err = VecDestroy(&sinfo.vector);PYLITH_CHECK_ERROR(err);
  err = DMCreateGlobalVector(sinfo.dm, &sinfo.vector);PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetName((PetscObject) sinfo.vector, _metadata.label.c_str());PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // createScatterWithBC

// ----------------------------------------------------------------------
// Create PETSc vector scatter for field. This is used to transfer
// information from the "global" PETSc vector view to the "local"
// PETSc section view. The PETSc vector includes constrained DOF. Use
// createScatter() if constrained DOF should be omitted from the PETSc
// vector.
void
pylith::topology::Field::createScatterWithBC(const Mesh& mesh,
					     const std::string& labelName,
					     PetscInt labelValue,
					     const char* context)
{ // createScatterWithBC
  PYLITH_METHOD_BEGIN;

  assert(context);
  PetscErrorCode err = 0;

  const bool createScatterOk = true;
  ScatterInfo& sinfo = _getScatter(context, createScatterOk);
  
  // Only create if scatter and scatterVec do not alreay exist.
  if (sinfo.dm) {
    assert(sinfo.vector);
    PYLITH_METHOD_END;
  } // if

  PetscDM dm = mesh.dmMesh();assert(dm);
  PetscSection section = NULL, newSection = NULL, gsection = NULL, subSection = NULL;
  PetscSF sf = NULL;
  PetscDMLabel subpointMap = NULL, subpointMapF = NULL;
  PetscInt dim, dimF, pStart, pEnd, qStart, qEnd, cEnd, cMax, vEnd, vMax;
  err = DMPlexGetHeightStratum(_dm, 0, NULL, &cEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetDepthStratum(_dm, 0, NULL, &vEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetHybridBounds(_dm, &cMax, NULL, NULL, &vMax);PYLITH_CHECK_ERROR(err);
  PetscInt excludeRanges[4] = {cMax, cEnd, vMax, vEnd};
  PetscInt numExcludes = (cMax >= 0 ? 1 : 0) + (vMax >= 0 ? 1 : 0);

  err = DMGetDefaultSection(_dm, &section);PYLITH_CHECK_ERROR(err);
  err = DMGetDimension(dm,  &dim);PYLITH_CHECK_ERROR(err);
  err = DMGetDimension(_dm, &dimF);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetChart(dm,  &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetChart(_dm, &qStart, &qEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetSubpointMap(dm,  &subpointMap);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetSubpointMap(_dm, &subpointMapF);PYLITH_CHECK_ERROR(err);
  if (((dim != dimF) || ((pEnd-pStart) < (qEnd-qStart))) && subpointMap && !subpointMapF) {
    const PetscInt *ind = NULL;
    PetscIS subpointIS = NULL;
    PetscInt n = 0, q = 0;

    err = PetscSectionGetChart(section, &qStart, &qEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexCreateSubpointIS(dm, &subpointIS);PYLITH_CHECK_ERROR(err);
    if (subpointIS) {
      err = ISGetLocalSize(subpointIS, &n);PYLITH_CHECK_ERROR(err);
      err = ISGetIndices(subpointIS, &ind);PYLITH_CHECK_ERROR(err);
    } // if
    err = PetscSectionCreate(mesh.comm(), &subSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetChart(subSection, pStart, pEnd);PYLITH_CHECK_ERROR(err);
    for(q = qStart; q < qEnd; ++q) {
      PetscInt dof, off, p;

      err = PetscSectionGetDof(section, q, &dof);PYLITH_CHECK_ERROR(err);
      if (dof) {
        err = PetscFindInt(q, n, ind, &p);PYLITH_CHECK_ERROR(err);
        if ((p >= pStart) && (p < pEnd)) {
          err = PetscSectionSetDof(subSection, p, dof);PYLITH_CHECK_ERROR(err);
          err = PetscSectionGetOffset(section, q, &off);PYLITH_CHECK_ERROR(err);
          err = PetscSectionSetOffset(subSection, p, off);PYLITH_CHECK_ERROR(err);
        } // if
      } // if
    } // for
    if (subpointIS) {
      err = ISRestoreIndices(subpointIS, &ind);PYLITH_CHECK_ERROR(err);
      err = ISDestroy(&subpointIS);PYLITH_CHECK_ERROR(err);
    } // if
    /* No need to setup section */
    section = subSection;
    /* There are no excludes for surface meshes */
    numExcludes = 0;
  } // if

  err = DMDestroy(&sinfo.dm);PYLITH_CHECK_ERROR(err);
  err = DMClone(dm, &sinfo.dm);PYLITH_CHECK_ERROR(err);
  err = PetscSectionClone(section, &newSection);PYLITH_CHECK_ERROR(err);
  err = DMSetDefaultSection(sinfo.dm, newSection);PYLITH_CHECK_ERROR(err);
  err = PetscSectionDestroy(&newSection);PYLITH_CHECK_ERROR(err);
  err = DMGetPointSF(sinfo.dm, &sf);PYLITH_CHECK_ERROR(err);
  if (labelName.empty()) {
    err = PetscSectionCreateGlobalSectionCensored(section, sf, PETSC_TRUE, numExcludes, excludeRanges, &gsection);PYLITH_CHECK_ERROR(err);
  } else {
    DMLabel label;

    err = DMPlexGetLabel(sinfo.dm, labelName.c_str(), &label);PYLITH_CHECK_ERROR(err);
    err = PetscSectionCreateGlobalSectionLabel(section, sf, PETSC_TRUE, label, labelValue, &gsection);PYLITH_CHECK_ERROR(err);
  } // if/else
  err = DMSetDefaultGlobalSection(sinfo.dm, gsection);PYLITH_CHECK_ERROR(err);
  err = PetscSectionDestroy(&gsection);PYLITH_CHECK_ERROR(err);
  err = VecDestroy(&sinfo.vector);PYLITH_CHECK_ERROR(err);
  err = DMCreateGlobalVector(sinfo.dm, &sinfo.vector);PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetName((PetscObject) sinfo.vector, _metadata.label.c_str());PYLITH_CHECK_ERROR(err);

  err = PetscSectionDestroy(&subSection);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // createScatterWithBC

// ----------------------------------------------------------------------
// Get PETSc vector associated with field.
PetscVec
pylith::topology::Field::vector(const char* context)
{ // vector
  PYLITH_METHOD_BEGIN;

  ScatterInfo& sinfo = _getScatter(context);

  PYLITH_METHOD_RETURN(sinfo.vector);
} // vector

// ----------------------------------------------------------------------
// Get PETSc vector associated with field.
const PetscVec
pylith::topology::Field::vector(const char* context) const
{ // vector
  PYLITH_METHOD_BEGIN;

  const ScatterInfo& sinfo = _getScatter(context);

  PYLITH_METHOD_RETURN(sinfo.vector);
} // vector

// ----------------------------------------------------------------------
// Scatter section information across processors to update the
//  PETSc vector view of the field.
void
pylith::topology::Field::scatterLocalToGlobal(const char* context) const
{ // scatterLocalToGlobal
  PYLITH_METHOD_BEGIN;

  assert(context);
  const ScatterInfo& sinfo = _getScatter(context);
  scatterLocalToGlobal(sinfo.vector, context);

  PYLITH_METHOD_END;
} // scatterLocalToGlobal

// ----------------------------------------------------------------------
// Scatter section information across processors to update the
//  PETSc vector view of the field.
void
pylith::topology::Field::scatterLocalToGlobal(const PetscVec vector,
						const char* context) const
{ // scatterLocalToGlobal
  PYLITH_METHOD_BEGIN;

  assert(vector);
  assert(context);
  const ScatterInfo& sinfo = _getScatter(context);
  PetscErrorCode err = 0;
  if (sinfo.dm) {
    err = DMLocalToGlobalBegin(sinfo.dm, _localVec, INSERT_VALUES, vector);PYLITH_CHECK_ERROR(err);
    err = DMLocalToGlobalEnd(sinfo.dm, _localVec, INSERT_VALUES, vector);PYLITH_CHECK_ERROR(err);
  } // if
  
  PYLITH_METHOD_END;
} // scatterLocalToGlobal

// ----------------------------------------------------------------------
// Scatter PETSc vector information across processors to update the
// section view of the field.
void
pylith::topology::Field::scatterGlobalToLocal(const char* context) const
{ // scatterGlobalToLocal
  PYLITH_METHOD_BEGIN;

  assert(context);

  const ScatterInfo& sinfo = _getScatter(context);
  scatterGlobalToLocal(sinfo.vector, context);

  PYLITH_METHOD_END;
} // scatterGlobalToLocal

// ----------------------------------------------------------------------
// Scatter PETSc vector information across processors to update the
// section view of the field.
void
pylith::topology::Field::scatterGlobalToLocal(const PetscVec vector,
					      const char* context) const
{ // scatterGlobalToLocal
  PYLITH_METHOD_BEGIN;

  assert(vector);
  assert(context);
  const ScatterInfo& sinfo = _getScatter(context);
  PetscErrorCode err = 0;

  if (sinfo.dm) {
    err = DMGlobalToLocalBegin(sinfo.dm, vector, INSERT_VALUES, _localVec);PYLITH_CHECK_ERROR(err);
    err = DMGlobalToLocalEnd(sinfo.dm, vector, INSERT_VALUES, _localVec);PYLITH_CHECK_ERROR(err);
  } // if

  PYLITH_METHOD_END;
} // scatterGlobalToLocal

// ----------------------------------------------------------------------
// Get scatter for given context.
pylith::topology::Field::ScatterInfo&
pylith::topology::Field::_getScatter(const char* context,
				     const bool createOk)
{ // _getScatter
  PYLITH_METHOD_BEGIN;

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
    err = DMDestroy(&sinfo.dm);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&sinfo.vector);PYLITH_CHECK_ERROR(err);

    _scatters.erase(context);
    isNewScatter = true;
  } // if

  if (isNewScatter && !createOk) {
    std::ostringstream msg;
    msg << "Scatter for context '" << context << "' does not exist for field '" << label() << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  ScatterInfo& sinfo = _scatters[context];
  if (isNewScatter) {
    sinfo.dm = 0;
    sinfo.vector = 0;
  } // if
  assert(_scatters.find(context) != _scatters.end());

  PYLITH_METHOD_RETURN(sinfo);
} // _getScatter

// ----------------------------------------------------------------------
// Get scatter for given context.
const pylith::topology::Field::ScatterInfo&
pylith::topology::Field::_getScatter(const char* context) const
{ // _getScatter
  PYLITH_METHOD_BEGIN;

  assert(context);

  const scatter_map_type::const_iterator s_iter = 
    _scatters.find(context);
  if (s_iter == _scatters.end()) {
    std::ostringstream msg;
    msg << "Scatter for context '" << context << "' does not exist for field '" << label() << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  PYLITH_METHOD_RETURN(s_iter->second);
} // _getScatter

// ----------------------------------------------------------------------
// Experimental
void
pylith::topology::Field::subfieldAdd(const char *name,
				     int numComponents,
				     const VectorFieldEnum fieldType,
				     const PylithScalar scale)
{ // subfieldAdd
  PYLITH_METHOD_BEGIN;

  assert(0 == _subfields.count(name));

  // Keep track of name/components until setup
  SubfieldInfo info;
  info.metadata.label = name;
  info.metadata.vectorFieldType = fieldType;
  info.metadata.scale = scale;
  info.metadata.dimsOkay = false;
  info.numComponents = numComponents;
  info.index = _subfields.size(); // Indices match order added.
  info.dm = NULL;
  _subfields[name] = info;
  
  PYLITH_METHOD_END;
} // subfieldAdd

// ----------------------------------------------------------------------
void
pylith::topology::Field::subfieldsSetup(void)
{ // subfieldsSetup
  PYLITH_METHOD_BEGIN;

  assert(_dm);

  // Setup section now that we know the total number of sub-fields and components.
  PetscSection section = NULL;
  PetscErrorCode err = DMGetDefaultSection(_dm, &section);PYLITH_CHECK_ERROR(err);assert(section);
  err = PetscSectionSetNumFields(section, _subfields.size());PYLITH_CHECK_ERROR(err);
  err = DMSetNumFields(_dm, _subfields.size());PYLITH_CHECK_ERROR(err);
  for (subfields_type::const_iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
    PetscObject fobj;
    const char* name = s_iter->first.c_str();
    const SubfieldInfo& info = s_iter->second;
    err = PetscSectionSetFieldName(section, info.index, name);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetFieldComponents(section, info.index, info.numComponents);PYLITH_CHECK_ERROR(err);
    err = DMGetField(_dm, info.index, &fobj);PYLITH_CHECK_ERROR(err);assert(section);
    err = PetscObjectSetName(fobj, name);PYLITH_CHECK_ERROR(err);
  } // for

  PYLITH_METHOD_END;
} // subfieldsSetup

// ----------------------------------------------------------------------
void
pylith::topology::Field::subfieldSetDof(const char *name,
					const DomainEnum domain,
					int fiberDim)
{ // subfieldSetDof
  PYLITH_METHOD_BEGIN;

  PetscInt pStart, pEnd;
  PetscErrorCode err;

  assert(_dm);
  switch(domain) {
  case VERTICES_FIELD:
    err = DMPlexGetDepthStratum(_dm, 0, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    break;
  case CELLS_FIELD:
    err = DMPlexGetHeightStratum(_dm, 0, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    break;
  case FACES_FIELD:
    err = DMPlexGetHeightStratum(_dm, 1, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    break;
  case POINTS_FIELD:
    err = DMPlexGetChart(_dm, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    break;
  default:
    std::ostringstream msg;
    msg << "Unknown value for DomainEnum: " << domain << "  in Field" << std::endl;
    throw std::logic_error(msg.str());
  } // switch
  PetscSection section = NULL;
  err = DMGetDefaultSection(_dm, &section);PYLITH_CHECK_ERROR(err);assert(section);
  const int iField = _subfields[name].index;
  for(PetscInt p = pStart; p < pEnd; ++p) {
    //err = PetscSectionAddDof(section, p, fiberDim);PYLITH_CHECK_ERROR(err); // Future use
    err = PetscSectionSetFieldDof(section, p, iField, fiberDim);PYLITH_CHECK_ERROR(err);
  } // for

  PYLITH_METHOD_END;
} // subfieldSetDof

// ----------------------------------------------------------------------
// Get metadata for subfield.
const pylith::topology::Field::SubfieldInfo&
pylith::topology::Field::subfieldInfo(const char* name) const
{ // subfieldInfo
  PYLITH_METHOD_BEGIN;

  subfields_type::const_iterator iter = _subfields.find(name);
  if (_subfields.end() == iter) {
    std::ostringstream msg;
    msg << "Could not find subfield '" << name << "' in field '" << label() << "'.\n"
	<< "Available subfields:";
    for (subfields_type::const_iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
      msg << " '" << s_iter->first << "'";
    } // for
    msg << std::endl;
    
    throw std::runtime_error(msg.str());
  } // if

  PYLITH_METHOD_RETURN(iter->second);
} // subfieldInfo


// ----------------------------------------------------------------------
// Copy subfield values to field.
void
pylith::topology::Field::copySubfield(const Field& field,
				      const char* name)
{ // copySubfield
  PYLITH_METHOD_BEGIN;

  // Check compatibility of sections
  const int srcSize = field.chartSize();
  const int dstSize = chartSize();
  if (dstSize < srcSize) {
    _extractSubfield(field, name);
  } // if
  assert(_localVec && field._localVec);

  const SubfieldInfo& subfieldInfo = const_cast<Field&>(field)._subfields[name];
  const int subfieldIndex = subfieldInfo.index;assert(subfieldIndex >= 0);

  _metadata = subfieldInfo.metadata;
  label(subfieldInfo.metadata.label.c_str()); // Use method to insure propagation to subsidiary objects

  PetscErrorCode err;
  const PetscSection& fieldSection = field.localSection();
  const PetscSection& subfieldSection = this->localSection();

  PetscInt pStart, pEnd;
  err = PetscSectionGetChart(subfieldSection, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

  // Copy values from field
  PylithScalar* subfieldArray = NULL;
  PylithScalar* fieldArray = NULL;
  err = VecGetArray(this->_localVec, &subfieldArray);PYLITH_CHECK_ERROR(err);
  err = VecGetArray(field._localVec, &fieldArray);PYLITH_CHECK_ERROR(err);
  for (PetscInt p = pStart; p < pEnd; ++p) {
    PetscInt fdof, foff, sdof, soff;
    
    err = PetscSectionGetDof(subfieldSection, p, &sdof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(subfieldSection, p, &soff);PYLITH_CHECK_ERROR(err);
    
    err = PetscSectionGetFieldDof(fieldSection, p, subfieldIndex, &fdof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetFieldOffset(fieldSection, p, subfieldIndex, &foff);PYLITH_CHECK_ERROR(err);
    
    assert(fdof == sdof);
    for (PetscInt d = 0; d < fdof; ++d) {
      subfieldArray[soff+d] = fieldArray[foff+d];
    } // for
  } // for
  err = VecRestoreArray(field._localVec, &fieldArray);PYLITH_CHECK_ERROR(err);
  err = VecRestoreArray(this->_localVec, &subfieldArray);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // copySubfield


// ----------------------------------------------------------------------
// Extract subfield.
void
pylith::topology::Field::_extractSubfield(const Field& field,
					  const char* name)
{ // _extractSubfield
  PYLITH_METHOD_BEGIN;

  clear();

  const SubfieldInfo& subfieldInfo = const_cast<Field&>(field)._subfields[name];
  const int subfieldIndex = subfieldInfo.index;assert(subfieldIndex >= 0);

  PetscErrorCode err;
  PetscIS subfieldIS = NULL;
  const int numSubfields = 1;
  int indicesSubfield[1];
  indicesSubfield[0] = subfieldIndex;
  err = DMDestroy(&_dm);PYLITH_CHECK_ERROR(err);
  if (subfieldInfo.dm) {
    err = DMClone(subfieldInfo.dm, &_dm);PYLITH_CHECK_ERROR(err);assert(_dm);
  } else {
    err = DMCreateSubDM(field.dmMesh(), numSubfields, indicesSubfield, &subfieldIS, &_dm);PYLITH_CHECK_ERROR(err);assert(_dm);
  } // if/else
  err = ISDestroy(&subfieldIS);PYLITH_CHECK_ERROR(err);
  
  err = DMCreateLocalVector(_dm, &_localVec);PYLITH_CHECK_ERROR(err);
  err = DMCreateGlobalVector(_dm, &_globalVec);PYLITH_CHECK_ERROR(err);

  // Setup section
  const PetscSection& fieldSection = field.localSection();
  PetscSection subfieldSection = NULL;
  err = DMGetDefaultSection(_dm, &subfieldSection);PYLITH_CHECK_ERROR(err);
  err = DMSetDefaultGlobalSection(_dm, NULL);PYLITH_CHECK_ERROR(err);

  PetscInt pStart = -1, pEnd = -1;
  err = PetscSectionGetChart(fieldSection, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
  err = PetscSectionSetChart(subfieldSection, pStart, pEnd);PYLITH_CHECK_ERROR(err);

  for (PetscInt p=pStart; p < pEnd; ++p) {
    PetscInt dof;
    err = PetscSectionGetFieldDof(fieldSection, p, subfieldIndex, &dof);PYLITH_CHECK_ERROR(err);
    if (dof > 0) {
      err = PetscSectionSetDof(subfieldSection, p, dof);PYLITH_CHECK_ERROR(err);

      err = PetscSectionGetFieldConstraintDof(fieldSection, p, subfieldIndex, &dof);PYLITH_CHECK_ERROR(err);
      err = PetscSectionSetConstraintDof(subfieldSection, p, dof);PYLITH_CHECK_ERROR(err);
    } // if
  } // for
  allocate();

  for (PetscInt p=pStart; p < pEnd; ++p) {
    PetscInt dof;
    const PetscInt* indices = NULL;
    err = PetscSectionGetConstraintDof(subfieldSection, p, &dof);PYLITH_CHECK_ERROR(err);
    if (dof > 0) {
      err = PetscSectionGetFieldConstraintIndices(fieldSection, p, subfieldIndex, &indices);PYLITH_CHECK_ERROR(err);
      err = PetscSectionSetConstraintIndices(subfieldSection, p, indices);PYLITH_CHECK_ERROR(err);
    } // if
  } // for

  PYLITH_METHOD_END;
} // _extractSubField

// End of file 
