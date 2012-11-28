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
template<typename mesh_type, typename section_type>
pylith::topology::Field<mesh_type, section_type>::Field(const mesh_type& mesh) :
  _mesh(mesh)
{ // constructor
  _metadata["default"].label = "unknown";
  _metadata["default"].vectorFieldType = OTHER;
  _metadata["default"].scale = 1.0;
  _metadata["default"].dimsOkay = false;
  if (mesh.dmMesh()) {
    DM             dm = mesh.dmMesh();
    Vec            coordVec;
    PetscSection   s;
    PetscErrorCode err;

    err = DMComplexClone(dm, &_dm);CHECK_PETSC_ERROR(err);
    err = DMGetCoordinatesLocal(dm, &coordVec);CHECK_PETSC_ERROR(err);
    if (coordVec) {
      DM           coordDM, newCoordDM;
      PetscSection coordSection, newCoordSection;

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
    _dm = PETSC_NULL;
  }
  _globalVec = PETSC_NULL;
  _localVec  = PETSC_NULL;
} // constructor

// ----------------------------------------------------------------------
// Constructor with mesh, section, and metadata.
template<typename mesh_type, typename section_type>
pylith::topology::Field<mesh_type, section_type>::Field(const mesh_type& mesh,
					  const ALE::Obj<section_type>& section,
					  const Metadata& metadata) :
  _mesh(mesh),
  _section(section)
{ // constructor
  assert(!section.isNull());
  _metadata["default"] = metadata;
  if (mesh.dmMesh()) {
    DM             dm = mesh.dmMesh(), coordDM, newCoordDM;
    PetscSection   coordSection, newCoordSection;
    Vec            coordVec;
    PetscSection   s;
    PetscErrorCode err;

    err = DMComplexClone(dm, &_dm);CHECK_PETSC_ERROR(err);
    err = DMGetCoordinatesLocal(dm, &coordVec);CHECK_PETSC_ERROR(err);
    if (coordVec) {
      err = DMGetCoordinateDM(dm, &coordDM);CHECK_PETSC_ERROR(err);
      err = DMGetCoordinateDM(_dm, &newCoordDM);CHECK_PETSC_ERROR(err);
      err = DMGetDefaultSection(coordDM, &coordSection);CHECK_PETSC_ERROR(err);
      err = PetscSectionClone(coordSection, &newCoordSection);CHECK_PETSC_ERROR(err);
      err = DMSetDefaultSection(newCoordDM, newCoordSection);CHECK_PETSC_ERROR(err);
      err = DMSetCoordinatesLocal(_dm, coordVec);CHECK_PETSC_ERROR(err);
    }
    this->dmSection(&s, &_localVec);
    err = DMSetDefaultSection(_dm, s);CHECK_PETSC_ERROR(err);
    err = DMCreateGlobalVector(_dm, &_globalVec);CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) _globalVec, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) _localVec,  _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  } else {
    _dm = PETSC_NULL;
    _globalVec = PETSC_NULL;
    _localVec  = PETSC_NULL;
  }
} // constructor

// ----------------------------------------------------------------------
// Constructor with mesh, DM, and metadata
template<typename mesh_type, typename section_type>
pylith::topology::Field<mesh_type, section_type>::Field(const mesh_type& mesh, DM dm, const Metadata& metadata) :
  _mesh(mesh),
  _dm(dm)
{ // constructor
  PetscErrorCode err;

  _metadata["default"] = metadata;
  err = DMCreateLocalVector(_dm, &_localVec);CHECK_PETSC_ERROR(err);
  err = DMCreateGlobalVector(_dm, &_globalVec);CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) _globalVec, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) _localVec,  _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
} // constructor

// ----------------------------------------------------------------------
// Constructor with field and subfields
template<typename mesh_type, typename section_type>
pylith::topology::Field<mesh_type, section_type>::Field(const Field& src,
                                                        const int fields[], int numFields) :
  _mesh(src._mesh),
  _section(PETSC_NULL)
{ // constructor
  DM             dm = mesh.dmMesh(), coordDM, newCoordDM;
  PetscSection   coordSection, newCoordSection;
  Vec            coordVec;
  PetscSection   s;
  PetscErrorCode err;

  _metadata["default"] = src._metadata["default"];
  err = DMGetDefaultSection(src._dm, &s);CHECK_PETSC_ERROR(err);
  for(PetscInt f = 0; f < numFields; ++f) {
    const char *name;

    err = PetscSectionGetFieldName(s, fields[f], &name);CHECK_PETSC_ERROR(err);
    _metadata[name] = src._metadata[name];
  }
  err = DMCreateSubDM(dm, numFields, fields, PETSC_NULL, &_dm);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dm, &coordVec);CHECK_PETSC_ERROR(err);
  if (coordVec) {
    err = DMGetCoordinateDM(dm, &coordDM);CHECK_PETSC_ERROR(err);
    err = DMGetCoordinateDM(_dm, &newCoordDM);CHECK_PETSC_ERROR(err);
    err = DMGetDefaultSection(coordDM, &coordSection);CHECK_PETSC_ERROR(err);
    err = PetscSectionClone(coordSection, &newCoordSection);CHECK_PETSC_ERROR(err);
    err = DMSetDefaultSection(newCoordDM, newCoordSection);CHECK_PETSC_ERROR(err);
    err = DMSetCoordinatesLocal(_dm, coordVec);CHECK_PETSC_ERROR(err);
  }
  _globalVec = PETSC_NULL;
  _localVec  = PETSC_NULL;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
template<typename mesh_type, typename section_type>
pylith::topology::Field<mesh_type, section_type>::~Field(void)
{ // destructor
  deallocate();
  PetscErrorCode err = DMDestroy(&_dm);CHECK_PETSC_ERROR(err);
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::deallocate(void)
{ // deallocate
  PetscErrorCode err = 0;
  
  const typename scatter_map_type::const_iterator scattersEnd = _scatters.end();
  for (typename scatter_map_type::iterator s_iter=_scatters.begin();
       s_iter != scattersEnd;
       ++s_iter) {

    err = DMDestroy(&s_iter->second.dm);CHECK_PETSC_ERROR(err);
    err = VecDestroy(&s_iter->second.vector);CHECK_PETSC_ERROR(err);

    if (s_iter->second.scatter) {err = VecScatterDestroy(&s_iter->second.scatter);CHECK_PETSC_ERROR(err);}
    err = VecDestroy(&s_iter->second.scatterVec);CHECK_PETSC_ERROR(err);
  } // for
  _scatters.clear();
  err = VecDestroy(&_globalVec);CHECK_PETSC_ERROR(err);
  err = VecDestroy(&_localVec);CHECK_PETSC_ERROR(err);
} // deallocate

// ----------------------------------------------------------------------
// Set label for field.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::label(const char* value)
{ // label
  _metadata["default"].label = value;
  if (!_section.isNull()) {
    _section->setName(value);
  } // if
  if (_localVec) {
    PetscErrorCode err = PetscObjectSetName((PetscObject) _localVec, value);CHECK_PETSC_ERROR(err);
  }
  if (_globalVec) {
    PetscErrorCode err = PetscObjectSetName((PetscObject) _globalVec, value);CHECK_PETSC_ERROR(err);
  }

  const typename scatter_map_type::const_iterator scattersEnd = _scatters.end();
  for (typename scatter_map_type::const_iterator s_iter=_scatters.begin();
       s_iter != scattersEnd;
       ++s_iter) {
    if (s_iter->second.vector) {
      PetscErrorCode err =
	PetscObjectSetName((PetscObject)s_iter->second.vector, value);
      CHECK_PETSC_ERROR(err);    
    } // if
  } // for
} // label

// ----------------------------------------------------------------------
// Get spatial dimension of domain.
template<typename mesh_type, typename section_type>
int
pylith::topology::Field<mesh_type, section_type>::spaceDim(void) const
{ // spaceDim
  const spatialdata::geocoords::CoordSys* cs = _mesh.coordsys();
  return (cs) ? cs->spaceDim() : 0;
} // spaceDim

// ----------------------------------------------------------------------
// Get the chart size.
template<typename mesh_type, typename section_type>
int
pylith::topology::Field<mesh_type, section_type>::chartSize(void) const
{ // chartSize
  if (_dm) {
    PetscSection   s;
    PetscInt       pStart, pEnd;
    PetscErrorCode err;

    err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetChart(s, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    return pEnd-pStart;
  }
  return _section.isNull() ? 0 : _section->getChart().size();
} // chartSize

// ----------------------------------------------------------------------
// Get the number of degrees of freedom.
template<typename mesh_type, typename section_type>
int
pylith::topology::Field<mesh_type, section_type>::sectionSize(void) const
{ // sectionSize
  return _section.isNull() ? 0 : _section->size();
} // sectionSize

// ----------------------------------------------------------------------
// Create seive section.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::newSection(void)
{ // newSection
  // Clear memory
  clear();

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");

  _section = new section_type(_mesh.comm(), _mesh.debug());  
  assert(!_section.isNull());
  _section->setName(_metadata["default"].label);

  logger.stagePop();
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion for a
// sequence of points.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::newSection(
				       const ALE::Obj<label_sequence>& points,
				       const int fiberDim)
{ // newSection
  typedef typename mesh_type::SieveMesh::point_type point_type;
  PetscErrorCode err;

  // Clear memory
  clear();

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");

  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _metadata["default"].label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  _section = new section_type(_mesh.comm(), _mesh.debug());
  assert(!_section.isNull());
  _section->setName(_metadata["default"].label);

  if (points->size() > 0) {
    const point_type pointMin = 
      *std::min_element(points->begin(), points->end());
    const point_type pointMax = 
      *std::max_element(points->begin(), points->end());
    _section->setChart(chart_type(pointMin, pointMax+1));
    _section->setFiberDimension(points, fiberDim);  

    if (_dm) {
      PetscSection s;
      err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetChart(s, pointMin, pointMax+1);CHECK_PETSC_ERROR(err);

      for(typename label_sequence::const_iterator p_iter = points->begin(); p_iter != points->end(); ++p_iter) {
        err = PetscSectionSetDof(s, *p_iter, fiberDim);CHECK_PETSC_ERROR(err);
      }
    }
  } else {// Create empty chart
    _section->setChart(chart_type(0, 0));
    if (_dm) {
      PetscSection s;
      err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetChart(s, 0, 0);CHECK_PETSC_ERROR(err);
    }
  }

  logger.stagePop();
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion for a list of
// points.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::newSection(const int_array& points,
					       const int fiberDim)
{ // newSection
  typedef typename mesh_type::SieveMesh::point_type point_type;
  PetscErrorCode err;

  // Clear memory
  clear();

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");
  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _metadata["default"].label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if
  
  _section = new section_type(_mesh.comm(), _mesh.debug());
  assert(!_section.isNull());
  _section->setName(_metadata["default"].label);

  const int npts = points.size();
  if (npts > 0) {
    const point_type pointMin = points.min();
    const point_type pointMax = points.max();
    _section->setChart(chart_type(pointMin, pointMax+1));
    for (int i=0; i < npts; ++i)
      _section->setFiberDimension(points[i], fiberDim);

    if (_dm) {
      PetscSection s;
      err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetChart(s, pointMin, pointMax+1);CHECK_PETSC_ERROR(err);

      for (int i=0; i < npts; ++i) {
        err = PetscSectionSetDof(s, points[i], fiberDim);CHECK_PETSC_ERROR(err);
      }
    }
  } else { // create empty chart
    _section->setChart(chart_type(0, 0));
    if (_dm) {
      PetscSection s;
      err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetChart(s, 0, 0);CHECK_PETSC_ERROR(err);
    }
  }

  logger.stagePop();
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::newSection(const DomainEnum domain,
					       const int fiberDim,
					       const int stratum)
{ // newSection
  // Changing this because cells/vertices are numbered differently in the new scheme
  if (_dm) {
    PetscSection   s;
    PetscInt       pStart, pEnd;
    PetscErrorCode err;

    switch(domain) {
    case VERTICES_FIELD:
      err = DMComplexGetDepthStratum(_dm, stratum, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
      break;
    case CELLS_FIELD:
      err = DMComplexGetHeightStratum(_dm, stratum, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
      break;
    case FACES_FIELD:
      err = DMComplexGetHeightStratum(_dm, stratum+1, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
      break;
    case POINTS_FIELD:
      err = DMComplexGetChart(_dm, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
      break;
    default:
      std::cerr << "Unknown value for DomainEnum: " << domain << std::endl;
      throw std::logic_error("Bad domain enum in Field.");
    }
    err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
    err = PetscSectionSetChart(s, pStart, pEnd);CHECK_PETSC_ERROR(err);

    for(PetscInt p = pStart; p < pEnd; ++p) {
      err = PetscSectionSetDof(s, p, fiberDim);CHECK_PETSC_ERROR(err);
    }
  } else {
  const ALE::Obj<SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  ALE::Obj<label_sequence> points;
  if (VERTICES_FIELD == domain)
    points = sieveMesh->depthStratum(stratum);
  else if (CELLS_FIELD == domain)
    points = sieveMesh->heightStratum(stratum);
  else if (FACES_FIELD == domain)
    points = sieveMesh->heightStratum(stratum+1);
  else {
    std::cerr << "Unknown value for DomainEnum: " << domain << std::endl;
    assert(0);
    throw std::logic_error("Bad domain enum in Field.");
  } // else

  newSection(points, fiberDim);
  }
} // newSection

// ----------------------------------------------------------------------
// Create section given chart.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::newSection(const Field& src,
					       const int fiberDim)
{ // newSection
  // Clear memory
  clear();

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");

  if (_section.isNull()) {
    logger.stagePop();
    newSection();
    logger.stagePush("Field");
  } // if
  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _metadata["default"].label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  const ALE::Obj<section_type>& srcSection = src._section;
  if (!srcSection.isNull()) {
    _section->setChart(srcSection->getChart());
    const chart_type& chart = _section->getChart();
    const typename chart_type::const_iterator chartBegin = chart.begin();
    const typename chart_type::const_iterator chartEnd = chart.end();

    for (typename chart_type::const_iterator c_iter = chartBegin;
	 c_iter != chartEnd;
	 ++c_iter) 
      if (srcSection->getFiberDimension(*c_iter) > 0)
	_section->setFiberDimension(*c_iter, fiberDim);

    if (_dm) {
      PetscSection   s;
      PetscInt       pStart = srcSection->getChart().min(), pEnd = srcSection->getChart().max();
      PetscErrorCode err;

      err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetChart(s, pStart, pEnd);CHECK_PETSC_ERROR(err);
      for(PetscInt p = pStart; p < pEnd; ++p) {
        err = PetscSectionSetDof(s, p, fiberDim);CHECK_PETSC_ERROR(err);
      }
    }
  } else if (src._dm) {
    PetscSection   srcs, s;
    PetscInt       pStart, pEnd;
    PetscErrorCode err;

    err = DMGetDefaultSection(src._dm, &srcs);CHECK_PETSC_ERROR(err);
    err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetChart(srcs, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    err = PetscSectionSetChart(s, pStart, pEnd);CHECK_PETSC_ERROR(err);
    for(PetscInt p = pStart; p < pEnd; ++p) {
      err = PetscSectionSetDof(s, p, fiberDim);CHECK_PETSC_ERROR(err);
    }
  } // if

  logger.stagePop();
} // newSection

// ----------------------------------------------------------------------
// Create section with same layout as another section.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::cloneSection(const Field& src)
{ // cloneSection
  std::string origLabel = _metadata["default"].label;

  // Clear memory
  clear();

  const ALE::Obj<section_type>& srcSection = src._section;
  if (!srcSection.isNull() && _section.isNull()) {
    newSection();
  } // if

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");

  _metadata["default"] = const_cast<Field&>(src)._metadata["default"];
  label(origLabel.c_str());

  if (!_section.isNull()) {
    if (!srcSection->sharedStorage()) {
      _section->setAtlas(srcSection->getAtlas());
      _section->allocateStorage();
      _section->setBC(srcSection->getBC());
      _section->copySpaces(srcSection);
    } else {
      _section->setChart(srcSection->getChart());
      const chart_type& chart = _section->getChart();
      const typename chart_type::const_iterator chartBegin = chart.begin();
      const typename chart_type::const_iterator chartEnd = chart.end();
      for (typename chart_type::const_iterator c_iter = chartBegin;
	   c_iter != chartEnd;
	   ++c_iter) {
        const int fiberDim = srcSection->getFiberDimension(*c_iter);
        if (fiberDim > 0)
          _section->setFiberDimension(*c_iter, fiberDim);
      } // for
      const ALE::Obj<SieveMesh>& sieveMesh = _mesh.sieveMesh();
      assert(!sieveMesh.isNull());
      sieveMesh->allocate(_section);
      _section->setBC(srcSection->getBC());
      _section->copySpaces(srcSection); // :BUG: NEED TO REBUILD SPACES 
    } // if/else
  } // if
  PetscSection   section = src.petscSection();
  PetscSection   newSection;
  PetscErrorCode err;

  if (_dm) {
    err = PetscSectionClone(section, &newSection);CHECK_PETSC_ERROR(err);
    err = DMSetDefaultSection(_dm, newSection);CHECK_PETSC_ERROR(err);
    err = DMCreateGlobalVector(_dm, &_globalVec);CHECK_PETSC_ERROR(err);
    err = DMCreateLocalVector(_dm, &_localVec);CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) _globalVec, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) _localVec,  _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  }
    
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
    }
      
	// Create scatter Vec
	if (s_iter->second.scatterVec) {
      int blockSize = 1;
      err = VecGetBlockSize(s_iter->second.scatterVec, &blockSize);CHECK_PETSC_ERROR(err);
      if (_section->sizeWithBC() > 0) {
        err = VecCreateSeqWithArray(PETSC_COMM_SELF,
                                    blockSize, _section->getStorageSize(),
                                    _section->restrictSpace(),
                                    &sinfo.scatterVec);
        CHECK_PETSC_ERROR(err);
      } else {
        err = VecCreateSeqWithArray(PETSC_COMM_SELF, 
                                    blockSize, 0, PETSC_NULL,
                                    &sinfo.scatterVec);
        CHECK_PETSC_ERROR(err);
      } // else
    }

	// Create vector using sizes from source section
	int vecLocalSize = 0;
	int vecGlobalSize = 0, vecGlobalSize2 = 0;
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
    }
    err = PetscObjectSetName((PetscObject)sinfo.vector, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
	
	_scatters[s_iter->first] = sinfo;
      } // for
    } // if
  logger.stagePop();
} // cloneSection

// ----------------------------------------------------------------------
// Get DMComplex section.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::dmSection(PetscSection *s, Vec *v) const {
  PetscSection   section;
  PetscInt       size;
  PetscInt       numNormalCells, numCohesiveCells, numNormalVertices, numShadowVertices, numLagrangeVertices;
  PetscErrorCode err;

  err = DMMeshConvertSection(_mesh.sieveMesh(), _section, &section);CHECK_PETSC_ERROR(err);
  _mesh.getPointTypeSizes(&numNormalCells, &numCohesiveCells, &numNormalVertices, &numShadowVertices, &numLagrangeVertices);
  if (numNormalCells+numCohesiveCells+numNormalVertices+numShadowVertices+numLagrangeVertices > 0) {
    PetscInt numFields, numComp, pMax, pStart, pEnd, qStart, qEnd;

    err = DMComplexGetChart(_mesh.dmMesh(), PETSC_NULL, &pMax);CHECK_PETSC_ERROR(err);
    err = PetscSectionCreate(_mesh.sieveMesh()->comm(), s);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetNumFields(section, &numFields);CHECK_PETSC_ERROR(err);
    if (numFields) {
      err = PetscSectionSetNumFields(*s, numFields);CHECK_PETSC_ERROR(err);
      for(PetscInt f = 0; f < numFields; ++f) {
        err = PetscSectionGetFieldComponents(section, f, &numComp);CHECK_PETSC_ERROR(err);
        err = PetscSectionSetFieldComponents(*s, f, numComp);CHECK_PETSC_ERROR(err);
      }
    }
    err = PetscSectionGetChart(section, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    if (pStart > 0) {
      qStart = pStart + numCohesiveCells;
    } else {
      qStart = pStart;
    }
    qEnd = PetscMin(pEnd + numCohesiveCells, pMax);
    err = PetscSectionSetChart(*s, qStart, qEnd);CHECK_PETSC_ERROR(err);
#if 0
    if (pEnd - pStart != (numNormalCells + numCohesiveCells + numNormalVertices + numShadowVertices + numLagrangeVertices)) {
      PetscPrintf(PETSC_COMM_SELF, "numCells %d != totCells %d\n", pEnd - pStart, numNormalCells + numCohesiveCells + numNormalVertices + numShadowVertices + numLagrangeVertices);
      CHECK_PETSC_ERROR(PETSC_ERR_ARG_SIZ);
    }
#endif
    /* The old-style point numbering: [normalCells, normalVertices, shadowVertices, lagrangeVertices, cohesiveCells]
       The new-style point numbering: [normalCells, cohesiveCells, normalVertices, shadowVertices, lagrangeVertices] */
    for(PetscInt p = pStart; p < numNormalCells; ++p) {
      PetscInt dof, fdof, cdof, cfdof, q = p;

      err = PetscSectionGetDof(section, p, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetDof(*s, q, dof);CHECK_PETSC_ERROR(err);
      for(PetscInt f = 0; f < numFields; ++f) {
        err = PetscSectionGetFieldDof(section, p, f, &fdof);CHECK_PETSC_ERROR(err);
        err = PetscSectionSetFieldDof(*s, q, f, fdof);CHECK_PETSC_ERROR(err);
      }
      err = PetscSectionGetConstraintDof(section, p, &cdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetConstraintDof(*s, q, cdof);CHECK_PETSC_ERROR(err);
      for(PetscInt f = 0; f < numFields; ++f) {
        err = PetscSectionGetFieldConstraintDof(section, p, f, &cfdof);CHECK_PETSC_ERROR(err);
        err = PetscSectionSetFieldConstraintDof(*s, q, f, cfdof);CHECK_PETSC_ERROR(err);
      }
    }
    for(PetscInt p = PetscMax(pStart, numNormalCells); p < PetscMin(pEnd, numNormalCells+numNormalVertices+numShadowVertices+numLagrangeVertices); ++p) {
      PetscInt dof, fdof, cdof, cfdof, q = p + numCohesiveCells;

      err = PetscSectionGetDof(section, p, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetDof(*s, q, dof);CHECK_PETSC_ERROR(err);
      for(PetscInt f = 0; f < numFields; ++f) {
        err = PetscSectionGetFieldDof(section, p, f, &fdof);CHECK_PETSC_ERROR(err);
        err = PetscSectionSetFieldDof(*s, q, f, fdof);CHECK_PETSC_ERROR(err);
      }
      err = PetscSectionGetConstraintDof(section, p, &cdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetConstraintDof(*s, q, cdof);CHECK_PETSC_ERROR(err);
      for(PetscInt f = 0; f < numFields; ++f) {
        err = PetscSectionGetFieldConstraintDof(section, p, f, &cfdof);CHECK_PETSC_ERROR(err);
        err = PetscSectionSetFieldConstraintDof(*s, q, f, cfdof);CHECK_PETSC_ERROR(err);
      }
    }
    for(PetscInt p = PetscMax(pStart, numNormalCells+numNormalVertices+numShadowVertices+numLagrangeVertices); p < pEnd; ++p) {
      PetscInt dof, fdof, cdof, cfdof, q = p - (numNormalVertices+numShadowVertices+numLagrangeVertices);

      err = PetscSectionGetDof(section, p, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetDof(*s, q, dof);CHECK_PETSC_ERROR(err);
      for(PetscInt f = 0; f < numFields; ++f) {
        err = PetscSectionGetFieldDof(section, p, f, &fdof);CHECK_PETSC_ERROR(err);
        err = PetscSectionSetFieldDof(*s, q, f, fdof);CHECK_PETSC_ERROR(err);
      }
      err = PetscSectionGetConstraintDof(section, p, &cdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetConstraintDof(*s, q, cdof);CHECK_PETSC_ERROR(err);
      for(PetscInt f = 0; f < numFields; ++f) {
        err = PetscSectionGetFieldConstraintDof(section, p, f, &cfdof);CHECK_PETSC_ERROR(err);
        err = PetscSectionSetFieldConstraintDof(*s, q, f, cfdof);CHECK_PETSC_ERROR(err);
      }
    }
    err = PetscSectionSetUp(*s);CHECK_PETSC_ERROR(err);
    for(PetscInt p = pStart; p < pStart+numNormalCells; ++p) {
      const PetscInt *cind, *cfind;
      PetscInt        q = p;

      err = PetscSectionGetConstraintIndices(section, p, &cind);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetConstraintIndices(*s, q, cind);CHECK_PETSC_ERROR(err);
      for(PetscInt f = 0; f < numFields; ++f) {
        err = PetscSectionGetFieldConstraintIndices(section, p, f, &cfind);CHECK_PETSC_ERROR(err);
        err = PetscSectionSetFieldConstraintIndices(*s, q, f, cfind);CHECK_PETSC_ERROR(err);
      }
    }
    for(PetscInt p = pStart+numNormalCells; p < pStart+numNormalCells+numNormalVertices+numShadowVertices+numLagrangeVertices; ++p) {
      const PetscInt *cind, *cfind;
      PetscInt        q = p + numCohesiveCells;

      err = PetscSectionGetConstraintIndices(section, p, &cind);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetConstraintIndices(*s, q, cind);CHECK_PETSC_ERROR(err);
      for(PetscInt f = 0; f < numFields; ++f) {
        err = PetscSectionGetFieldConstraintIndices(section, p, f, &cfind);CHECK_PETSC_ERROR(err);
        err = PetscSectionSetFieldConstraintIndices(*s, q, f, cfind);CHECK_PETSC_ERROR(err);
      }
    }
    for(PetscInt p = pStart+numNormalCells+numNormalVertices+numShadowVertices+numLagrangeVertices; p < pEnd; ++p) {
      const PetscInt *cind, *cfind;
      PetscInt        q = p - (numNormalVertices+numShadowVertices+numLagrangeVertices);

      err = PetscSectionGetConstraintIndices(section, p, &cind);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetConstraintIndices(*s, q, cind);CHECK_PETSC_ERROR(err);
      for(PetscInt f = 0; f < numFields; ++f) {
        err = PetscSectionGetFieldConstraintIndices(section, p, f, &cfind);CHECK_PETSC_ERROR(err);
        err = PetscSectionSetFieldConstraintIndices(*s, q, f, cfind);CHECK_PETSC_ERROR(err);
      }
    }
    err = PetscSectionDestroy(&section);CHECK_PETSC_ERROR(err);
  } else {
    *s = section;
  }
  err = PetscSectionGetStorageSize(*s, &size);CHECK_PETSC_ERROR(err);
  err = VecCreateMPIWithArray(_section->comm(), 1, size, PETSC_DETERMINE, _section->restrictSpace(), v);CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) *v, _section->getName().c_str());CHECK_PETSC_ERROR(err);
}

// ----------------------------------------------------------------------
// Clear variables associated with section.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::clear(void)
{ // clear
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");

  deallocate();
  if (!_section.isNull())
    _section->clear();

  _metadata["default"].scale = 1.0;
  _metadata["default"].vectorFieldType = OTHER;
  _metadata["default"].dimsOkay = false;

  logger.stagePop();
} // clear

// ----------------------------------------------------------------------
// Allocate Sieve section.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::allocate(void)
{ // allocate
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");
  PetscSection   s = PETSC_NULL;
  PetscErrorCode err;

  if (_dm) {err = DMGetDefaultSection(_dm, &s);CHECK_PETSC_ERROR(err);}
  assert(!_section.isNull() || s);

  if (!_section.isNull()) {
    const ALE::Obj<SieveMesh>& sieveMesh = _mesh.sieveMesh();
    assert(!sieveMesh.isNull());
    sieveMesh->allocate(_section);
  }
  if (s) {
    err = PetscSectionSetUp(s);CHECK_PETSC_ERROR(err);
    err = DMCreateGlobalVector(_dm, &_globalVec);CHECK_PETSC_ERROR(err);
    err = DMCreateLocalVector(_dm, &_localVec);CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) _globalVec, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject) _localVec,  _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  }

  logger.stagePop();
} // allocate

// ----------------------------------------------------------------------
// Zero section values (excluding constrained DOF).
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::zero(void)
{ // zero
  if (!_section.isNull())
    _section->zero(); // Does not zero BC.
  if (_localVec) {
    PetscSection   section;
    PetscInt       pStart, pEnd, maxDof = 0;
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
        err = DMComplexVecSetClosure(_dm, section, _localVec, p, &values[0], INSERT_VALUES);CHECK_PETSC_ERROR(err);
      } // if
    } // for
  }
} // zero

// ----------------------------------------------------------------------
// Zero section values (including constrained DOF).
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::zeroAll(void)
{ // zeroAll
  if (!_section.isNull()) {
    const chart_type& chart = _section->getChart();
    const typename chart_type::const_iterator chartBegin = chart.begin();
    const typename chart_type::const_iterator chartEnd = chart.end();

    // Assume fiber dimension is uniform
    const int fiberDim = (chart.size() > 0) ? 
      _section->getFiberDimension(*chartBegin) : 0;
    scalar_array values(fiberDim);
    values *= 0.0;

    for (typename chart_type::const_iterator c_iter = chartBegin;
	 c_iter != chartEnd;
	 ++c_iter) {
      if (_section->getFiberDimension(*c_iter) > 0) {
	assert(fiberDim == _section->getFiberDimension(*c_iter));
	_section->updatePointAll(*c_iter, &values[0]);
      } // if
    } // for
  } // if
  if (_localVec) {
    PetscErrorCode err = VecSet(_localVec, 0.0);CHECK_PETSC_ERROR(err);
  }
} // zeroAll

// ----------------------------------------------------------------------
// Complete section by assembling across processors.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::complete(void)
{ // complete
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Completion");

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  if (!_section.isNull())
    ALE::Completion::completeSectionAdd(sieveMesh->getSendOverlap(),
					sieveMesh->getRecvOverlap(), 
					_section, _section);

  if (_dm) {
    // Not sure if DMLocalToLocal() would work
    PetscErrorCode err;

    err = VecSet(_globalVec, 0.0);CHECK_PETSC_ERROR(err);
    err = DMLocalToGlobalBegin(_dm, _localVec, ADD_VALUES, _globalVec);CHECK_PETSC_ERROR(err);
    err = DMLocalToGlobalEnd(_dm, _localVec, ADD_VALUES, _globalVec);CHECK_PETSC_ERROR(err);
    err = DMGlobalToLocalBegin(_dm, _globalVec, INSERT_VALUES, _localVec);CHECK_PETSC_ERROR(err);
    err = DMGlobalToLocalEnd(_dm, _globalVec, INSERT_VALUES, _localVec);CHECK_PETSC_ERROR(err);
  }
  logger.stagePop();
} // complete

// ----------------------------------------------------------------------
// Copy field values and metadata.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::copy(const Field& field)
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

  if (!_section.isNull() && !field._section.isNull()) {
    // Copy values from field
    const chart_type& chart = _section->getChart();
    const typename chart_type::const_iterator chartBegin = chart.begin();
    const typename chart_type::const_iterator chartEnd = chart.end();

    for (typename chart_type::const_iterator c_iter = chartBegin;
	 c_iter != chartEnd;
	 ++c_iter) {
      assert(field._section->getFiberDimension(*c_iter) ==
	     _section->getFiberDimension(*c_iter));
      if (_section->getFiberDimension(*c_iter) > 0)
	_section->updatePointAll(*c_iter, 
				 field._section->restrictPoint(*c_iter));
    } // for
  } // if

  if (_localVec && field._localVec) {
    PetscErrorCode err = VecCopy(field._localVec, _localVec);CHECK_PETSC_ERROR(err);
  }

  label(const_cast<Field&>(field)._metadata["default"].label.c_str()); // Update label
  _metadata["default"].scale = const_cast<Field&>(field)._metadata["default"].scale;
} // copy

// ----------------------------------------------------------------------
// Copy field values.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::copy(const ALE::Obj<section_type>& osection)
{ // copy
  // Check compatibility of sections
  const int srcSize = osection->getChart().size();
  const int dstSize = chartSize();
  if (srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot copy values from Sieve section "
	<< _metadata["default"].label << "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    size: " << srcSize << "\n"
	<< "  Destination section:\n"
	<< "    space dim: " << spaceDim() << "\n"
	<< "    vector field type: " << _metadata["default"].vectorFieldType << "\n"
	<< "    scale: " << _metadata["default"].scale << "\n"
	<< "    size: " << dstSize;
    throw std::runtime_error(msg.str());
  } // if
  assert( (_section.isNull() && osection.isNull()) ||
	  (!_section.isNull() && !osection.isNull()) );

  if (!_section.isNull()) {
    // Copy values from field
    const chart_type& chart = _section->getChart();
    const typename chart_type::const_iterator chartEnd = chart.end();

    for (typename chart_type::const_iterator c_iter = chart.begin();
	 c_iter != chartEnd;
	 ++c_iter) {
      assert(osection->getFiberDimension(*c_iter) ==
	     _section->getFiberDimension(*c_iter));
      if (osection->getFiberDimension(*c_iter))
	_section->updatePoint(*c_iter, osection->restrictPoint(*c_iter));
    } // for
  } // if
} // copy

template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::copy(PetscSection osection, PetscInt field, PetscInt component, Vec ovec)
{ // copy
  PetscSection   section;
  PetscScalar   *array, *oarray;
  PetscInt       numFields, numComp, pStart, pEnd, qStart, qEnd;
  PetscErrorCode err;

  assert(_dm);
  err = DMGetDefaultSection(_dm, &section);CHECK_PETSC_ERROR(err);
  assert(section);assert(_localVec);
  assert(osection);assert(ovec);
  err = PetscSectionGetNumFields(section, &numFields);CHECK_PETSC_ERROR(err);
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
  if (field >= 0) {
    err = PetscSectionGetFieldComponents(section, field, &numComp);CHECK_PETSC_ERROR(err);
    if (component >= numComp) {
      std::ostringstream msg;
      msg << "Invalid field component "<<component<<" should be in [0, "<<numComp<<")";
      throw std::runtime_error(msg.str());
    }
  }
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
      assert(!(odof%numComp));
      odof  = odof/numComp;
      ooff += odof*component;
    } else {
      err = PetscSectionGetDof(osection, p, &odof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(osection, p, &ooff);CHECK_PETSC_ERROR(err);
    }
    assert(odof == dof);
    if (!odof) continue;
    for(PetscInt d = 0; d < dof; ++d) {
      array[off+d] = oarray[ooff+d];
    }
  } // for
  err = VecRestoreArray(_localVec, &array);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(ovec, &oarray);CHECK_PETSC_ERROR(err);
} // copy

// ----------------------------------------------------------------------
// Add two fields, storing the result in one of the fields.
template<typename mesh_type, typename section_type>
pylith::topology::Field<mesh_type, section_type>&
pylith::topology::Field<mesh_type, section_type>::operator+=(const Field& field)
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
  assert( (_section.isNull() && field._section.isNull()) ||
	  (!_section.isNull() && !field._section.isNull()) );

  if (!_section.isNull()) {
    // Add values from field
    const chart_type& chart = _section->getChart();
    const typename chart_type::const_iterator chartBegin = chart.begin();
    const typename chart_type::const_iterator chartEnd = chart.end();

    // Assume fiber dimension is uniform
    const int fiberDim = (chart.size() > 0) ? 
      _section->getFiberDimension(*chartBegin) : 0;
    scalar_array values(fiberDim);

    for (typename chart_type::const_iterator c_iter = chartBegin;
	 c_iter != chartEnd;
	 ++c_iter) {
      if (field._section->getFiberDimension(*c_iter) > 0) {
	assert(fiberDim == field._section->getFiberDimension(*c_iter));
	assert(fiberDim == _section->getFiberDimension(*c_iter));
	field._section->restrictPoint(*c_iter, &values[0], values.size());
	_section->updatePointAllAdd(*c_iter, &values[0]);
      } // if
    } // for
  } // if

  if (_localVec && field._localVec) {
    PetscErrorCode err = VecAXPY(_localVec, 1.0, field._localVec);CHECK_PETSC_ERROR(err);
  }

  return *this;
} // operator+=

// ----------------------------------------------------------------------
// Dimensionalize field.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::dimensionalize(void) const
{ // dimensionalize
  if (!const_cast<Field*>(this)->_metadata["default"].dimsOkay) {
    std::ostringstream msg;
    msg << "Cannot dimensionalize field '" << const_cast<Field*>(this)->_metadata["default"].label
	<< "' because the flag "
	<< "has been set to keep field nondimensional.";
    throw std::runtime_error(msg.str());
  } // if

  spatialdata::units::Nondimensional normalizer;

  if (!_section.isNull()) {
    const chart_type& chart = _section->getChart();
    const typename chart_type::const_iterator chartEnd = chart.end();

    // Assume fiber dimension is uniform
    const int fiberDim = (chart.size() > 0) ? 
      _section->getFiberDimension(*chart.begin()) : 0;
    scalar_array values(fiberDim);

    for (typename chart_type::const_iterator c_iter = chart.begin();
	 c_iter != chartEnd;
	 ++c_iter) 
      if (_section->getFiberDimension(*c_iter) > 0) {
	assert(fiberDim == _section->getFiberDimension(*c_iter));
      
	_section->restrictPoint(*c_iter, &values[0], values.size());
	normalizer.dimensionalize(&values[0], values.size(), const_cast<Field*>(this)->_metadata["default"].scale);
	_section->updatePointAll(*c_iter, &values[0]);
      } // if
  } // if

  if (_localVec) {
    PetscSection   section;
    PetscScalar   *array;
    PetscInt       pStart, pEnd;
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
  }
} // dimensionalize

// ----------------------------------------------------------------------
// Print field to standard out.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::view(const char* label) const
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
  if (!_section.isNull())
    _section->view(label);
  if (_dm) {
    PetscSection   section;
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
template<typename mesh_type, typename section_type>
template<typename scatter_mesh_type>
void
pylith::topology::Field<mesh_type, section_type>::createScatter(const scatter_mesh_type& mesh,
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

  if (!_section.isNull()) {
    assert(!mesh.sieveMesh().isNull());
    // Get global order (create if necessary).
    const std::string& orderLabel = _section->getName();
    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());
    const ALE::Obj<typename mesh_type::SieveMesh::order_type>& order = 
      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, orderLabel, _section);
    assert(!order.isNull());

    // Create scatter
    err = DMMeshCreateGlobalScatter(sieveMesh, _section, order, false, &sinfo.scatter);
    CHECK_PETSC_ERROR(err);
  
    // Create scatterVec
    const int blockSize = 1;
    if (_section->sizeWithBC() > 0) {
      err = VecCreateSeqWithArray(PETSC_COMM_SELF,
                                  blockSize, _section->getStorageSize(),
                                  _section->restrictSpace(),
                                  &sinfo.scatterVec);CHECK_PETSC_ERROR(err);
    } else {
      err = VecCreateSeqWithArray(PETSC_COMM_SELF, 
                                  blockSize, 0, PETSC_NULL,
                                  &sinfo.scatterVec);CHECK_PETSC_ERROR(err);
    } // else

#if 0
    // Create vector
    err = VecCreate(_mesh.comm(), &sinfo.vector);CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject)sinfo.vector, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
    err = VecSetSizes(sinfo.vector, order->getLocalSize(), order->getGlobalSize());CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(sinfo.vector, blockSize);CHECK_PETSC_ERROR(err);
    err = VecSetFromOptions(sinfo.vector);CHECK_PETSC_ERROR(err);
#endif
  }
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
template<typename mesh_type, typename section_type>
template<typename scatter_mesh_type>
void
pylith::topology::Field<mesh_type, section_type>::createScatter(
      const scatter_mesh_type& mesh,
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

  if (!_section.isNull()); {
    assert(!mesh.sieveMesh().isNull());
    // Get global order (create if necessary).
    const std::string& orderLabel = 
      (strlen(context) > 0) ?
      _section->getName() + std::string("_") + std::string(context) :
      _section->getName();
    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());
    const ALE::Obj<typename mesh_type::SieveMesh::order_type>& order = 
      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, orderLabel,
                                              numbering->getChart().begin(),
                                              numbering->getChart().end(),
                                              _section);
    assert(!order.isNull());

    // Create scatter
    err = DMMeshCreateGlobalScatter(sieveMesh, _section, order, false, &sinfo.scatter);CHECK_PETSC_ERROR(err);

    // Create scatterVec
    const int blockSize = 1;
    if (_section->sizeWithBC() > 0) {
      err = VecCreateSeqWithArray(PETSC_COMM_SELF,
                                  blockSize, _section->getStorageSize(),
                                  _section->restrictSpace(),
                                  &sinfo.scatterVec);CHECK_PETSC_ERROR(err);
    } else {
      err = VecCreateSeqWithArray(PETSC_COMM_SELF, 
                                  blockSize, 0, PETSC_NULL,
                                  &sinfo.scatterVec);CHECK_PETSC_ERROR(err);
    } // else

#if 0
    // Create vector
    err = VecCreate(mesh.comm(), &sinfo.vector);CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject)sinfo.vector, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
    err = VecSetSizes(sinfo.vector, order->getLocalSize(), order->getGlobalSize());CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(sinfo.vector, blockSize);CHECK_PETSC_ERROR(err);
    err = VecSetFromOptions(sinfo.vector); CHECK_PETSC_ERROR(err);  
#endif
  }
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
template<typename mesh_type, typename section_type>
template<typename scatter_mesh_type>
void
pylith::topology::Field<mesh_type, section_type>::createScatterWithBC(
        const scatter_mesh_type& mesh,
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

  if (!_section.isNull()) {
    assert(!mesh.sieveMesh().isNull());

    // Get global order (create if necessary).
    const std::string& orderLabel = _section->getName();
    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());
    const ALE::Obj<typename mesh_type::SieveMesh::order_type>& order = 
      sieveMesh->getFactory()->getGlobalOrderWithBC(sieveMesh, orderLabel, _section);
    assert(!order.isNull());

    // Create scatter
    err = DMMeshCreateGlobalScatter(sieveMesh, _section, order, true, &sinfo.scatter); 
    CHECK_PETSC_ERROR(err);
  
    // Create scatterVec
    const int blockSize = _getFiberDim();
    if (_section->sizeWithBC() > 0) {
      err = VecCreateSeqWithArray(PETSC_COMM_SELF,
                                  blockSize, _section->getStorageSize(),
                                  _section->restrictSpace(),
                                  &sinfo.scatterVec);CHECK_PETSC_ERROR(err);
    } else {
      err = VecCreateSeqWithArray(PETSC_COMM_SELF, 
                                  blockSize, 0, PETSC_NULL,
                                  &sinfo.scatterVec);CHECK_PETSC_ERROR(err);
    } // else
#if 0
    // Create vector
    err = VecCreate(mesh.comm(), &sinfo.vector);CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject)sinfo.vector, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
    err = VecSetSizes(sinfo.vector, order->getLocalSize(), order->getGlobalSize());CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(sinfo.vector, blockSize);CHECK_PETSC_ERROR(err);
    err = VecSetFromOptions(sinfo.vector);CHECK_PETSC_ERROR(err);
#endif
  }

  PetscSection section, newSection, gsection;
  PetscSF      sf;

  err = DMComplexClone(_dm, &sinfo.dm);CHECK_PETSC_ERROR(err);
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
template<typename mesh_type, typename section_type>
template<typename scatter_mesh_type>
void
pylith::topology::Field<mesh_type, section_type>::createScatterWithBC(
       const scatter_mesh_type& mesh,
       const typename ALE::Obj<typename SieveMesh::numbering_type> numbering,
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

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("GlobalOrder");

  if (!_section.isNull() && !numbering.isNull()) {
    // Get global order (create if necessary).
    const std::string& orderLabel = 
      (strlen(context) > 0) ?
      _section->getName() + std::string("_") + std::string(context) :
      _section->getName();
    const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());
    const ALE::Obj<typename mesh_type::SieveMesh::order_type>& order = 
      sieveMesh->getFactory()->getGlobalOrderWithBC(sieveMesh, orderLabel,
                                                    numbering->getChart().begin(),
                                                    numbering->getChart().end(),
                                                    _section);
    assert(!order.isNull());
    //order->view("GLOBAL ORDER"); // DEBUG

    // Create scatter
    err = DMMeshCreateGlobalScatter(sieveMesh, _section, order, true, &sinfo.scatter); 
    CHECK_PETSC_ERROR(err);

    // Create scatterVec
    const int blockSize = _getFiberDim();
    if (_section->sizeWithBC() > 0) {
      err = VecCreateSeqWithArray(PETSC_COMM_SELF,
                                  blockSize, _section->getStorageSize(),
                                  _section->restrictSpace(),
                                  &sinfo.scatterVec);CHECK_PETSC_ERROR(err);
    } else {
      err = VecCreateSeqWithArray(PETSC_COMM_SELF, 
                                  blockSize, 0, PETSC_NULL,
                                  &sinfo.scatterVec);CHECK_PETSC_ERROR(err);
    } // else

#if 0
    // Create vector
    err = VecCreate(mesh.comm(), &sinfo.vector);CHECK_PETSC_ERROR(err);
    err = PetscObjectSetName((PetscObject)sinfo.vector, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
    err = VecSetSizes(sinfo.vector,order->getLocalSize(), order->getGlobalSize());CHECK_PETSC_ERROR(err);
    err = VecSetBlockSize(sinfo.vector, blockSize);CHECK_PETSC_ERROR(err);
    err = VecSetFromOptions(sinfo.vector); CHECK_PETSC_ERROR(err);  
#endif
  }

  PetscSection section, newSection, gsection;
  PetscSF      sf;
  PetscInt     cEnd, cMax, vEnd, vMax;
  err = DMComplexGetHeightStratum(_dm, 0, PETSC_NULL, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMComplexGetDepthStratum(_dm, 0, PETSC_NULL, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMComplexGetVTKBounds(_dm, &cMax, &vMax);CHECK_PETSC_ERROR(err);
  PetscInt     excludeRanges[4] = {cMax, cEnd, vMax, vEnd};
  PetscInt     numExcludes      = (cMax >= 0 ? 1 : 0) + (vMax >= 0 ? 1 : 0);

  err = DMComplexClone(_dm, &sinfo.dm);CHECK_PETSC_ERROR(err);
  err = DMGetDefaultSection(_dm, &section);CHECK_PETSC_ERROR(err);
  err = PetscSectionClone(section, &newSection);CHECK_PETSC_ERROR(err);
  err = DMSetDefaultSection(sinfo.dm, newSection);CHECK_PETSC_ERROR(err);
  err = DMGetPointSF(sinfo.dm, &sf);CHECK_PETSC_ERROR(err);
  err = PetscSectionCreateGlobalSectionCensored(section, sf, PETSC_TRUE, numExcludes, excludeRanges, &gsection);CHECK_PETSC_ERROR(err);
  err = DMSetDefaultGlobalSection(sinfo.dm, gsection);CHECK_PETSC_ERROR(err);
  err = DMCreateGlobalVector(sinfo.dm, &sinfo.vector);CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject) sinfo.vector, _metadata["default"].label.c_str());CHECK_PETSC_ERROR(err);
  PetscInt localSize, globalSize;

  err = PetscSectionGetStorageSize(section, &localSize);CHECK_PETSC_ERROR(err);
  err = VecGetSize(sinfo.vector, &globalSize);CHECK_PETSC_ERROR(err);
  /* assert(order->getLocalSize()  == localSize); This does not work because the local vector includes the lagrange cell variables */
  /* assert(order->getGlobalSize() == globalSize); */

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
  
  logger.stagePop();
} // createScatterWithBC

// ----------------------------------------------------------------------
// Get PETSc vector associated with field.
template<typename mesh_type, typename section_type>
PetscVec
pylith::topology::Field<mesh_type, section_type>::vector(const char* context)
{ // vector
  std::ostringstream msg;

  ScatterInfo& sinfo = _getScatter(context);
  return sinfo.vector;
} // vector

// ----------------------------------------------------------------------
// Get PETSc vector associated with field.
template<typename mesh_type, typename section_type>
const PetscVec
pylith::topology::Field<mesh_type, section_type>::vector(const char* context) const
{ // vector
  std::ostringstream msg;

  const ScatterInfo& sinfo = _getScatter(context);
  return sinfo.vector;
} // vector

// ----------------------------------------------------------------------
// Scatter section information across processors to update the
//  PETSc vector view of the field.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::scatterSectionToVector(const char* context) const
{ // scatterSectionToVector
  assert(context);

  const ScatterInfo& sinfo = _getScatter(context);
  scatterSectionToVector(sinfo.vector, context);
} // scatterSectionToVector

// ----------------------------------------------------------------------
// Scatter section information across processors to update the
//  PETSc vector view of the field.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::scatterSectionToVector(const PetscVec vector,
									 const char* context) const
{ // scatterSectionToVector
  assert(vector);
  assert(context);
  const ScatterInfo& sinfo = _getScatter(context);
  PetscErrorCode     err   = 0;
#if 0
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
  }
} // scatterSectionToVector

// ----------------------------------------------------------------------
// Scatter PETSc vector information across processors to update the
// section view of the field.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::scatterVectorToSection(const char* context) const
{ // scatterVectorToSection
  assert(context);

  const ScatterInfo& sinfo = _getScatter(context);
  scatterVectorToSection(sinfo.vector, context);
} // scatterVectorToSection

// ----------------------------------------------------------------------
// Scatter PETSc vector information across processors to update the
// section view of the field.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::scatterVectorToSection(const PetscVec vector,
									 const char* context) const
{ // scatterVectorToSection
  assert(vector);
  assert(context);
  const ScatterInfo& sinfo = _getScatter(context);
  PetscErrorCode     err   = 0;

#if 0
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
  }
} // scatterVectorToSection

// ----------------------------------------------------------------------
// Setup split field with all one space per spatial dimension.
template<typename mesh_type, typename section_type>
void 
pylith::topology::Field<mesh_type, section_type>::splitDefault(void)
{ // splitDefault
  assert(!_section.isNull());
  const int spaceDim = _mesh.dimension();
  for (int iDim=0; iDim < spaceDim; ++iDim)
    _section->addSpace(); // displacements

  const chart_type& chart = _section->getChart();

  const typename chart_type::const_iterator chartBegin = chart.begin();
  const typename chart_type::const_iterator chartEnd = chart.end();
  for (int fibration=0; fibration < spaceDim; ++fibration)
    for (typename chart_type::const_iterator c_iter = chart.begin();
        c_iter != chartEnd;
        ++c_iter) {
      if (_section->getFiberDimension(*c_iter) > 0) {
	assert(spaceDim == _section->getFiberDimension(*c_iter));
	_section->setFiberDimension(*c_iter, 1, fibration);
      } // if
    } // for
} // splitDefault

// ----------------------------------------------------------------------
// Get fiber dimension associated with section (only works if fiber
// dimension is uniform).
template<typename mesh_type, typename section_type>
int
pylith::topology::Field<mesh_type, section_type>::_getFiberDim(void)
{ // _getFiberDim
  
  int fiberDimLocal = (this->chartSize() > 0) ? _section->getFiberDimension(*_section->getChart().begin()) : 0;
  int fiberDim = 0;
  MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPI_INT, MPI_MAX,
		_mesh.comm());

  return fiberDim;
} // _getFiberDim

// ----------------------------------------------------------------------
// Get scatter for given context.
template<typename mesh_type, typename section_type>
typename pylith::topology::Field<mesh_type, section_type>::ScatterInfo&
pylith::topology::Field<mesh_type, section_type>::_getScatter(const char* context,
							      const bool createOk)
{ // _getScatter
  assert(context);

  bool isNewScatter = _scatters.find(context) == _scatters.end();

  // Synchronize creation of scatter (empty sections may have
  // leftover, reusable scatters that need to be cleared out).
  int numNewScatterLocal = (isNewScatter) ? 1 : 0;
  int numNewScatter = 0;
  MPI_Allreduce(&numNewScatterLocal, &numNewScatter, 1, MPI_INT, MPI_MAX,
		_mesh.comm());
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
template<typename mesh_type, typename section_type>
const typename pylith::topology::Field<mesh_type, section_type>::ScatterInfo&
pylith::topology::Field<mesh_type, section_type>::_getScatter(const char* context) const
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

// Experimental
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::addField(const char *name, int numComponents)
{
  // Keep track of name/components until setup
  _tmpFields[name] = numComponents;
  _metadata[name]  = _metadata["default"];
}

template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::setupFields()
{
  assert(_dm);
  // Keep track of name/components until setup
  PetscSection   section;
  PetscInt       f = 0;
  PetscErrorCode err = DMGetDefaultSection(_dm, &section);CHECK_PETSC_ERROR(err);
  assert(section);
  err = PetscSectionSetNumFields(section, _tmpFields.size());CHECK_PETSC_ERROR(err);
  for(std::map<std::string, int>::const_iterator f_iter = _tmpFields.begin(); f_iter != _tmpFields.end(); ++f_iter, ++f) {
    err = PetscSectionSetFieldName(section, f, f_iter->first.c_str());CHECK_PETSC_ERROR(err);
    err = PetscSectionSetFieldComponents(section, f, f_iter->second);CHECK_PETSC_ERROR(err);
  }
  _tmpFields.clear();
#if 0
  // Right now, we assume that the section covers the entire chart
  PetscInt pStart, pEnd;

  err = DMComplexGetChart(_dm, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
  err = PetscSectionSetChart(section, pStart, pEnd);CHECK_PETSC_ERROR(err);
#endif
}

template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::updateDof(const char *name, const DomainEnum domain, int fiberDim)
{
  PetscSection   section;
  PetscInt       pStart, pEnd, f = 0;
  PetscErrorCode err;

  assert(_dm);
  switch(domain) {
  case VERTICES_FIELD:
    err = DMComplexGetDepthStratum(_dm, 0, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    break;
  case CELLS_FIELD:
    err = DMComplexGetHeightStratum(_dm, 0, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    break;
  case FACES_FIELD:
    err = DMComplexGetHeightStratum(_dm, 1, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    break;
  case POINTS_FIELD:
    err = DMComplexGetChart(_dm, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    break;
  default:
    std::cerr << "Unknown value for DomainEnum: " << domain << std::endl;
    throw std::logic_error("Bad domain enum in Field.");
  }
  err = DMGetDefaultSection(_dm, &section);CHECK_PETSC_ERROR(err);
  assert(section);
  for(map_type::const_iterator f_iter = _metadata.begin(); f_iter != _metadata.end(); ++f_iter) {
    if (f_iter->first == name) break;
    if (f_iter->first == "default") continue;
    ++f;
  }
  assert(f < _metadata.size());
  for(PetscInt p = pStart; p < pEnd; ++p) {
    //err = PetscSectionAddDof(section, p, fiberDim);CHECK_PETSC_ERROR(err);
    err = PetscSectionSetFieldDof(section, p, f, fiberDim);CHECK_PETSC_ERROR(err);
  }
}

// End of file 
