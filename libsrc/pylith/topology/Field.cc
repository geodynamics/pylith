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

#include "pylith/utils/array.hh" // USES double_array

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
  _metadata.label = "unknown";
  _metadata.vectorFieldType = OTHER;
  _metadata.scale = 1.0;
  _metadata.dimsOkay = false;
} // constructor

// ----------------------------------------------------------------------
// Constructor with mesh, section, and metadata.
template<typename mesh_type, typename section_type>
pylith::topology::Field<mesh_type, section_type>::Field(const mesh_type& mesh,
					  const ALE::Obj<section_type>& section,
					  const Metadata& metadata) :
  _metadata(metadata),
  _mesh(mesh),
  _section(section)
{ // constructor
  assert(!section.isNull());
} // constructor

// ----------------------------------------------------------------------
// Destructor.
template<typename mesh_type, typename section_type>
pylith::topology::Field<mesh_type, section_type>::~Field(void)
{ // destructor
  deallocate();
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
    std::cout << "SCATTER DESTROY, label: " << label() << ", context: " << s_iter->first << std::endl;

    err = VecDestroy(&s_iter->second.vector);CHECK_PETSC_ERROR(err);
    err = VecScatterDestroy(&s_iter->second.scatter);CHECK_PETSC_ERROR(err);
    err = VecDestroy(&s_iter->second.scatterVec);CHECK_PETSC_ERROR(err);
  } // for
  _scatters.clear();
} // deallocate

// ----------------------------------------------------------------------
// Set label for field.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::label(const char* value)
{ // label
  _metadata.label = value;
  if (!_section.isNull()) {
    _section->setName(value);
  } // if

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
  _section->setName(_metadata.label);

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

  // Clear memory
  clear();

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");

  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _metadata.label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  _section = new section_type(_mesh.comm(), _mesh.debug());
  assert(!_section.isNull());
  _section->setName(_metadata.label);

  if (points->size() > 0) {
    const point_type pointMin = 
      *std::min_element(points->begin(), points->end());
    const point_type pointMax = 
      *std::max_element(points->begin(), points->end());
    _section->setChart(chart_type(pointMin, pointMax+1));
    _section->setFiberDimension(points, fiberDim);  
  } else // Create empty chart
    _section->setChart(chart_type(0, 0));

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

  // Clear memory
  clear();

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");
  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _metadata.label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if
  
  _section = new section_type(_mesh.comm(), _mesh.debug());
  assert(!_section.isNull());
  _section->setName(_metadata.label);

  const int npts = points.size();
  if (npts > 0) {
    const point_type pointMin = points.min();
    const point_type pointMax = points.max();
    _section->setChart(chart_type(pointMin, pointMax+1));
    for (int i=0; i < npts; ++i)
      _section->setFiberDimension(points[i], fiberDim);
  } else  // create empty chart
    _section->setChart(chart_type(0, 0));

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
  const ALE::Obj<SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  ALE::Obj<label_sequence> points;
  if (VERTICES_FIELD == domain)
    points = sieveMesh->depthStratum(stratum);
  else if (CELLS_FIELD == domain)
    points = sieveMesh->heightStratum(stratum);
  else {
    std::cerr << "Unknown value for DomainEnum: " << domain << std::endl;
    assert(0);
    throw std::logic_error("Bad domain enum in Field.");
  } // else

  newSection(points, fiberDim);
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
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _metadata.label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  const ALE::Obj<section_type>& srcSection = src.section();
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
  } // if

  logger.stagePop();
} // newSection

// ----------------------------------------------------------------------
// Create section with same layout as another section.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::cloneSection(const Field& src)
{ // cloneSection
  std::string origLabel = _metadata.label;

  // Clear memory
  clear();

  const ALE::Obj<section_type>& srcSection = src.section();
  if (!srcSection.isNull() && _section.isNull()) {
    newSection();
  } // if

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");

  _metadata = src._metadata;
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
    
    // Reuse scatters in clone
    PetscErrorCode err = 0;
    if (src._scatters.size() > 0) {
      const typename scatter_map_type::const_iterator scattersEnd = src._scatters.end();
      for (typename scatter_map_type::const_iterator s_iter=src._scatters.begin();
	   s_iter != scattersEnd;
	   ++s_iter) {
	ScatterInfo sinfo;
	sinfo.vector = 0;
	sinfo.scatterVec = 0;

    std::cout << "SCATTER COPY, label: " << label() << ", context: " << s_iter->first << std::endl;

	// Copy scatter
	sinfo.scatter = s_iter->second.scatter;
	err = PetscObjectReference((PetscObject) sinfo.scatter);
	CHECK_PETSC_ERROR(err);
      
	// Create scatter Vec
	if (_section->sizeWithBC() > 0) {
	  err = VecCreateSeqWithArray(PETSC_COMM_SELF,
				      _section->getStorageSize(),
				      _section->restrictSpace(),
				      &sinfo.scatterVec);
	  CHECK_PETSC_ERROR(err);
	} else {
	  err = VecCreateSeqWithArray(PETSC_COMM_SELF, 0, 0,
				      &sinfo.scatterVec);
	  CHECK_PETSC_ERROR(err);
	} // else

	// Create vector using sizes from source section
	int vecLocalSize = 0;
	int vecGlobalSize = 0;
	err = VecGetLocalSize(s_iter->second.vector, &vecLocalSize); 
	CHECK_PETSC_ERROR(err);
	err = VecGetSize(s_iter->second.vector, &vecGlobalSize); CHECK_PETSC_ERROR(err);

	err = VecCreate(_mesh.comm(), &sinfo.vector); CHECK_PETSC_ERROR(err);
	err = PetscObjectSetName((PetscObject)sinfo.vector, _metadata.label.c_str());
	CHECK_PETSC_ERROR(err);
	err = VecSetSizes(sinfo.vector, vecLocalSize, vecGlobalSize); 
	CHECK_PETSC_ERROR(err);
	err = VecSetFromOptions(sinfo.vector); CHECK_PETSC_ERROR(err);  
	
	_scatters[s_iter->first] = sinfo;
      } // for
    } // if
  } // if
  logger.stagePop();
} // cloneSection

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

  _metadata.scale = 1.0;
  _metadata.vectorFieldType = OTHER;
  _metadata.dimsOkay = false;

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

  assert(!_section.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  sieveMesh->allocate(_section);

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
    double_array values(fiberDim);
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
      field._metadata.vectorFieldType != _metadata.vectorFieldType ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot copy values from section '" << field._metadata.label 
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
  assert( (_section.isNull() && field._section.isNull()) ||
	  (!_section.isNull() && !field._section.isNull()) );

  if (!_section.isNull()) {
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

  label(field._metadata.label.c_str()); // Update label
  _metadata.scale = field._metadata.scale;
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
	<< _metadata.label << "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    size: " << srcSize << "\n"
	<< "  Destination section:\n"
	<< "    space dim: " << spaceDim() << "\n"
	<< "    vector field type: " << _metadata.vectorFieldType << "\n"
	<< "    scale: " << _metadata.scale << "\n"
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
    double_array values(fiberDim);

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

  return *this;
} // operator+=

// ----------------------------------------------------------------------
// Dimensionalize field.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::dimensionalize(void) const
{ // dimensionalize
  if (!_metadata.dimsOkay) {
    std::ostringstream msg;
    msg << "Cannot dimensionalize field '" << _metadata.label
	<< "' because the flag "
	<< "has been set to keep field nondimensional.";
    throw std::runtime_error(msg.str());
  } // if

  if (!_section.isNull()) {
    const chart_type& chart = _section->getChart();
    const typename chart_type::const_iterator chartEnd = chart.end();

    // Assume fiber dimension is uniform
    const int fiberDim = (chart.size() > 0) ? 
      _section->getFiberDimension(*chart.begin()) : 0;
    double_array values(fiberDim);

    spatialdata::units::Nondimensional normalizer;

    for (typename chart_type::const_iterator c_iter = chart.begin();
	 c_iter != chartEnd;
	 ++c_iter) 
      if (_section->getFiberDimension(*c_iter) > 0) {
	assert(fiberDim == _section->getFiberDimension(*c_iter));
      
	_section->restrictPoint(*c_iter, &values[0], values.size());
	normalizer.dimensionalize(&values[0], values.size(), _metadata.scale);
	_section->updatePointAll(*c_iter, &values[0]);
      } // if
  } // if
} // dimensionalize

// ----------------------------------------------------------------------
// Print field to standard out.
template<typename mesh_type, typename section_type>
void
pylith::topology::Field<mesh_type, section_type>::view(const char* label) const
{ // view
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
    default :
      std::cerr << "Unknown vector field value '" << _metadata.vectorFieldType
		<< "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad vector field type in Field.");
    } // switch

  std::cout << "Viewing field '" << _metadata.label << "' "<< label << ".\n"
	    << "  vector field type: " << vecFieldString << "\n"
	    << "  scale: " << _metadata.scale << "\n"
	    << "  dimensionalize flag: " << _metadata.dimsOkay << std::endl;
  if (!_section.isNull())
    _section->view(label);
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
  assert(!_section.isNull());
  assert(!mesh.sieveMesh().isNull());

  PetscErrorCode err = 0;
  const bool createScatterOk = true;
  ScatterInfo& sinfo = _getScatter(context, createScatterOk);
  if (sinfo.scatter) {
    assert(sinfo.scatterVec);
    assert(sinfo.vector);
    return;
  } // if

    std::cout << "SCATTER CREATE, label: " << label() << ", context: " << context << std::endl;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("GlobalOrder");

  // Get global order (create if necessary).
  const std::string& orderLabel = _section->getName();
  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh =
    mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<typename mesh_type::SieveMesh::order_type>& order = 
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, orderLabel,
					    _section);
  assert(!order.isNull());

  // Create scatter
  err = DMMeshCreateGlobalScatter(sieveMesh, _section, order, false, 
				  &sinfo.scatter);
  CHECK_PETSC_ERROR(err);
  
  // Create scatterVec
  if (_section->sizeWithBC() > 0)  {
    err = VecCreateSeqWithArray(PETSC_COMM_SELF, _section->getStorageSize(),
				_section->restrictSpace(),
				&sinfo.scatterVec); CHECK_PETSC_ERROR(err);
  } else {
    err = VecCreateSeqWithArray(PETSC_COMM_SELF, 0, 0,
				&sinfo.scatterVec); CHECK_PETSC_ERROR(err);
  } // if/else

  // Create vector
  err = VecCreate(_mesh.comm(), &sinfo.vector);
  CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject)sinfo.vector,
			   _metadata.label.c_str());
  CHECK_PETSC_ERROR(err);
  err = VecSetSizes(sinfo.vector,
		    order->getLocalSize(), order->getGlobalSize());
  CHECK_PETSC_ERROR(err);
  err = VecSetFromOptions(sinfo.vector); CHECK_PETSC_ERROR(err);  
  
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
  assert(!_section.isNull());
  assert(!mesh.sieveMesh().isNull());

  PetscErrorCode err = 0;
  const bool createScatterOk = true;
  ScatterInfo& sinfo = _getScatter(context, createScatterOk);
  
  // Only create if scatter and scatterVec do not alreay exist.
  if (sinfo.scatter) {
    assert(sinfo.scatterVec);
    assert(sinfo.vector);
    return;
  } // if

  std::cout << "SCATTER CREATE, label: " << label() << ", context: " << context << std::endl;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("GlobalOrder");

  // Get global order (create if necessary).
  const std::string& orderLabel = 
    (strlen(context) > 0) ?
    _section->getName() + std::string("_") + std::string(context) :
    _section->getName();
  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh =
    mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<typename mesh_type::SieveMesh::order_type>& order = 
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, orderLabel,
                                            numbering->getChart().begin(),
                                            numbering->getChart().end(),
                                            _section);
  assert(!order.isNull());

  // Create scatter
  err = DMMeshCreateGlobalScatter(sieveMesh, _section, order, false, 
				  &sinfo.scatter); 
  CHECK_PETSC_ERROR(err);

  // Create scatterVec
  if (_section->sizeWithBC() > 0)  {
    err = VecCreateSeqWithArray(PETSC_COMM_SELF, _section->getStorageSize(),
				_section->restrictSpace(),
				&sinfo.scatterVec); CHECK_PETSC_ERROR(err);
  } else {
    err = VecCreateSeqWithArray(PETSC_COMM_SELF, 0, 0,
				&sinfo.scatterVec); CHECK_PETSC_ERROR(err);
  } // if/else

  // Create vector
  err = VecCreate(mesh.comm(), &sinfo.vector);
  CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject)sinfo.vector,
			   _metadata.label.c_str());
  CHECK_PETSC_ERROR(err);
  err = VecSetSizes(sinfo.vector,
		    order->getLocalSize(), order->getGlobalSize());
  CHECK_PETSC_ERROR(err);
  err = VecSetFromOptions(sinfo.vector); CHECK_PETSC_ERROR(err);  

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
  assert(!_section.isNull());
  assert(!mesh.sieveMesh().isNull());


  PetscErrorCode err = 0;
  const bool createScatterOk = true;
  ScatterInfo& sinfo = _getScatter(context, createScatterOk);
  if (sinfo.scatter) {
    assert(sinfo.scatterVec);
    assert(sinfo.vector);
    return;
  } // if

    std::cout << "SCATTER CREATE, label: " << label() << ", context: " << context << std::endl;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("GlobalOrder");

  // Get global order (create if necessary).
  const std::string& orderLabel = _section->getName();
  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh =
    mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<typename mesh_type::SieveMesh::order_type>& order = 
    sieveMesh->getFactory()->getGlobalOrderWithBC(sieveMesh, orderLabel,
						  _section);
  assert(!order.isNull());

  // Create scatter
  err = DMMeshCreateGlobalScatter(sieveMesh, _section, order, true, 
				  &sinfo.scatter); 
  CHECK_PETSC_ERROR(err);
  
  // Create scatterVec
  if (_section->sizeWithBC() > 0)  {
    err = VecCreateSeqWithArray(PETSC_COMM_SELF, _section->getStorageSize(),
				_section->restrictSpace(),
				&sinfo.scatterVec); CHECK_PETSC_ERROR(err);
  } else {
    err = VecCreateSeqWithArray(PETSC_COMM_SELF, 0, 0,
				&sinfo.scatterVec); CHECK_PETSC_ERROR(err);
  } // if/else
  
  // Create vector
   err = VecCreate(mesh.comm(), &sinfo.vector);
  CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject)sinfo.vector,
			   _metadata.label.c_str());
  CHECK_PETSC_ERROR(err);
  err = VecSetSizes(sinfo.vector,
		    order->getLocalSize(), order->getGlobalSize());
  CHECK_PETSC_ERROR(err);
  err = VecSetFromOptions(sinfo.vector); CHECK_PETSC_ERROR(err);  
  
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
  assert(!numbering.isNull());
  assert(context);
  assert(!_section.isNull());

  PetscErrorCode err = 0;
  const bool createScatterOk = true;
  ScatterInfo& sinfo = _getScatter(context, createScatterOk);
  
  // Only create if scatter and scatterVec do not alreay exist.
  if (sinfo.scatter) {
    assert(sinfo.scatterVec);
    assert(sinfo.vector);
    return;
  } // if

    std::cout << "SCATTER CREATE, label: " << label() << ", context: " << context << std::endl;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("GlobalOrder");

  // Get global order (create if necessary).
  const std::string& orderLabel = 
    (strlen(context) > 0) ?
    _section->getName() + std::string("_") + std::string(context) :
    _section->getName();
  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh =
    mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<typename mesh_type::SieveMesh::order_type>& order = 
    sieveMesh->getFactory()->getGlobalOrderWithBC(sieveMesh, orderLabel,
                                                  numbering->getChart().begin(),
                                                  numbering->getChart().end(),
                                                  _section);
  assert(!order.isNull());
  //order->view("GLOBAL ORDER"); // DEBUG

  // Create scatter
  err = DMMeshCreateGlobalScatter(sieveMesh, _section, order, true, 
				  &sinfo.scatter); 
  CHECK_PETSC_ERROR(err);

  // Create scatterVec
  if (_section->sizeWithBC() > 0)  {
    err = VecCreateSeqWithArray(PETSC_COMM_SELF, _section->getStorageSize(),
				_section->restrictSpace(),
				&sinfo.scatterVec); CHECK_PETSC_ERROR(err);
  } else {
    err = VecCreateSeqWithArray(PETSC_COMM_SELF, 0, 0,
				&sinfo.scatterVec); CHECK_PETSC_ERROR(err);
  } // if/else

  // Create vector
  err = VecCreate(mesh.comm(), &sinfo.vector);
  CHECK_PETSC_ERROR(err);
  err = PetscObjectSetName((PetscObject)sinfo.vector,
			   _metadata.label.c_str());
  CHECK_PETSC_ERROR(err);
  err = VecSetSizes(sinfo.vector,
		    order->getLocalSize(), order->getGlobalSize());
  CHECK_PETSC_ERROR(err);
  err = VecSetFromOptions(sinfo.vector); CHECK_PETSC_ERROR(err);  

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
  ScatterInfo& sinfo = _getScatter(context);
  return sinfo.vector;
} // vector

// ----------------------------------------------------------------------
// Get PETSc vector associated with field.
template<typename mesh_type, typename section_type>
const PetscVec
pylith::topology::Field<mesh_type, section_type>::vector(const char* context) const
{ // vector
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
  assert(!_section.isNull());

  PetscErrorCode err = 0;
  const ScatterInfo& sinfo = _getScatter(context);
  err = VecScatterBegin(sinfo.scatter, sinfo.scatterVec, vector,
			INSERT_VALUES, SCATTER_FORWARD); CHECK_PETSC_ERROR(err);
  err = VecScatterEnd(sinfo.scatter, sinfo.scatterVec, vector,
		      INSERT_VALUES, SCATTER_FORWARD); CHECK_PETSC_ERROR(err);
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
  assert(!_section.isNull());

  PetscErrorCode err = 0;
  const ScatterInfo& sinfo = _getScatter(context);
  err = VecScatterBegin(sinfo.scatter, vector, sinfo.scatterVec,
			INSERT_VALUES, SCATTER_REVERSE); CHECK_PETSC_ERROR(err);
  err = VecScatterEnd(sinfo.scatter, vector, sinfo.scatterVec,
		      INSERT_VALUES, SCATTER_REVERSE); CHECK_PETSC_ERROR(err);
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


// End of file 
