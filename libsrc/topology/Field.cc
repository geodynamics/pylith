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
// Copyright (c) 2010 University of California, Davis
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
template<typename mesh_type>
pylith::topology::Field<mesh_type>::Field(const mesh_type& mesh) :
  _scale(1.0),
  _label("unknown"),
  _mesh(mesh),
  _vector(0),
  _scatter(0),
  _scatterVec(0),
  _vecFieldType(OTHER),
  _dimensionsOkay(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
template<typename mesh_type>
pylith::topology::Field<mesh_type>::~Field(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::deallocate(void)
{ // deallocate
  PetscErrorCode err = 0;
  if (0 != _vector) {
    err = VecDestroy(&_vector); _vector = 0;
    CHECK_PETSC_ERROR(err);
  } // if

  if (0 != _scatter) {
    err = VecScatterDestroy(&_scatter); _scatter = 0;
    CHECK_PETSC_ERROR(err);
  } // if

  if (0 != _scatterVec) {
    err = VecDestroy(&_scatterVec); _scatterVec = 0;
    CHECK_PETSC_ERROR(err);
  } // if
} // deallocate

// ----------------------------------------------------------------------
// Get spatial dimension of domain.
template<typename mesh_type>
int
pylith::topology::Field<mesh_type>::spaceDim(void) const
{ // spaceDim
  const spatialdata::geocoords::CoordSys* cs = _mesh.coordsys();
  return (0 != cs) ? cs->spaceDim() : 0;
} // spaceDim

// ----------------------------------------------------------------------
// Get the chart size.
template<typename mesh_type>
int
pylith::topology::Field<mesh_type>::chartSize(void) const
{ // chartSize
  return _section.isNull() ? 0 : _section->getChart().size();
} // chartSize

// ----------------------------------------------------------------------
// Get the number of degrees of freedom.
template<typename mesh_type>
int
pylith::topology::Field<mesh_type>::sectionSize(void) const
{ // sectionSize
  return _section.isNull() ? 0 : _section->size();
} // sectionSize

// ----------------------------------------------------------------------
// Create seive section.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(void)
{ // newSection
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //std::cout << "Making Field " << _label << " empty section" << std::endl;
  logger.stagePush("Field");
  _section = new RealSection(_mesh.comm(), _mesh.debug());  
  logger.stagePop();
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion for a
// sequence of points.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(
				       const ALE::Obj<label_sequence>& points,
				       const int fiberDim)
{ // newSection
  typedef typename mesh_type::SieveMesh::point_type point_type;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //std::cout << "Making Field " << _label << " section type 1" << std::endl;
  logger.stagePush("Field");
  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  _section = new RealSection(_mesh.comm(), _mesh.debug());

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
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(const int_array& points,
					       const int fiberDim)
{ // newSection
  typedef typename mesh_type::SieveMesh::point_type point_type;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //std::cout << "Making Field " << _label << " section type 1b" << std::endl;
  logger.stagePush("Field");
  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if
  
  _section = new RealSection(_mesh.comm(), _mesh.debug());

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
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(const DomainEnum domain,
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
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(const Field& src,
					       const int fiberDim)
{ // newSection
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //std::cout << "Making Field " << _label << " section type 2" << std::endl;
  logger.stagePush("Field");
  if (_section.isNull()) {
    logger.stagePop();
    newSection();
    logger.stagePush("Field");
  } // if
  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for field '" << _label
	<< "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  const ALE::Obj<RealSection>& srcSection = src.section();
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

  //std::cout << "Done making Field " << _label << " section type 2" << std::endl;
  logger.stagePop();
} // newSection

// ----------------------------------------------------------------------
// Create section with same layout as another section.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::cloneSection(const Field& src)
{ // cloneSection
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //std::cout << "Making Field " << _label << " section type 3" << std::endl;
  logger.stagePush("Field");

  deallocate();
  _vecFieldType = src._vecFieldType;
  _scale = src._scale;
  _dimensionsOkay = false;

  const ALE::Obj<RealSection>& srcSection = src.section();
  if (!srcSection.isNull() && _section.isNull()) {
    logger.stagePop();
    newSection();
    logger.stagePush("Field");
  }

  if (!_section.isNull()) {
    _section->setAtlas(srcSection->getAtlas());
    _section->allocateStorage();
    _section->setBC(srcSection->getBC());
    _section->copySpaces(srcSection);

    PetscErrorCode err = 0;
    if (0 != src._scatter) {
      _scatter = src._scatter;
      err = PetscObjectReference((PetscObject) _scatter);
      CHECK_PETSC_ERROR(err);

      if (_section->sizeWithBC() > 0) {
	err = VecCreateSeqWithArray(PETSC_COMM_SELF,
				    _section->getStorageSize(),
				    _section->restrictSpace(),
				    &_scatterVec);
	CHECK_PETSC_ERROR(err);
      } else {
	err = VecCreateSeqWithArray(PETSC_COMM_SELF, 0, 0,
				    &_scatterVec);
	CHECK_PETSC_ERROR(err);
      } // else
    } // if
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
  if (!_section.isNull())
    _section->clear();

  _scale = 1.0;
  _vecFieldType = OTHER;
  _dimensionsOkay = false;

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

  assert(!_section.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  sieveMesh->allocate(_section);

  logger.stagePop();
} // allocate

// ----------------------------------------------------------------------
// Zero section values (excluding constrained DOF).
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::zero(void)
{ // zero
  if (!_section.isNull())
    _section->zero(); // Does not zero BC.
} // zero

// ----------------------------------------------------------------------
// Zero section values (including constrained DOF).
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::zeroAll(void)
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
      if (0 != _section->getFiberDimension(*c_iter)) {
	assert(fiberDim == _section->getFiberDimension(*c_iter));
	_section->updatePointAll(*c_iter, &values[0]);
      } // if
    } // for
  } // if
} // zeroAll

// ----------------------------------------------------------------------
// Complete section by assembling across processors.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::complete(void)
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
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::copy(const Field& field)
{ // copy
  // Check compatibility of sections
  const int srcSize = field.chartSize();
  const int dstSize = chartSize();
  if (field.spaceDim() != spaceDim() ||
      field._vecFieldType != _vecFieldType ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot copy values from section '" << field._label 
	<< "' to section '" << _label << "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    space dim: " << field.spaceDim() << "\n"
	<< "    vector field type: " << field._vecFieldType << "\n"
	<< "    scale: " << field._scale << "\n"
	<< "    size: " << srcSize << "\n"
	<< "  Destination section:\n"
	<< "    space dim: " << spaceDim() << "\n"
	<< "    vector field type: " << _vecFieldType << "\n"
	<< "    scale: " << _scale << "\n"
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
      _section->updatePointAll(*c_iter, field._section->restrictPoint(*c_iter));
    } // for
  } // if

  _label = field._label;
  _scale = field._scale;
} // copy

// ----------------------------------------------------------------------
// Copy field values.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::copy(const ALE::Obj<typename mesh_type::RealSection>& osection)
{ // copy
  // Check compatibility of sections
  const int srcSize = osection->getChart().size();
  const int dstSize = chartSize();
  if (srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot copy values from Sieve section "
	<< _label << "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    size: " << srcSize << "\n"
	<< "  Destination section:\n"
	<< "    space dim: " << spaceDim() << "\n"
	<< "    vector field type: " << _vecFieldType << "\n"
	<< "    scale: " << _scale << "\n"
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
      _section->updatePoint(*c_iter, osection->restrictPoint(*c_iter));
    } // for
  } // if
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
      field._vecFieldType != _vecFieldType ||
      field._scale != _scale ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot add values from section '" << field._label 
	<< "' to section '" << _label << "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    space dim: " << field.spaceDim() << "\n"
	<< "    vector field type: " << field._vecFieldType << "\n"
	<< "    scale: " << field._scale << "\n"
	<< "    size: " << srcSize << "\n"
	<< "  Destination section:\n"
	<< "    space dim: " << spaceDim() << "\n"
	<< "    vector field type: " << _vecFieldType << "\n"
	<< "    scale: " << _scale << "\n"
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
      if (0 != field._section->getFiberDimension(*c_iter)) {
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
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::dimensionalize(void) const
{ // dimensionalize
  if (!_dimensionsOkay) {
    std::ostringstream msg;
    msg << "Cannot dimensionalize field '" << _label << "' because the flag "
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
      if (0 != _section->getFiberDimension(*c_iter)) {
	assert(fiberDim == _section->getFiberDimension(*c_iter));
      
	_section->restrictPoint(*c_iter, &values[0], values.size());
	normalizer.dimensionalize(&values[0], values.size(), _scale);
	_section->updatePointAll(*c_iter, &values[0]);
      } // if
  } // if
} // dimensionalize

// ----------------------------------------------------------------------
// Print field to standard out.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::view(const char* label) const
{ // view
  std::string vecFieldString;
  switch(_vecFieldType)
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
      std::cerr << "Unknown vector field value '" << _vecFieldType
		<< "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad vector field type in Field.");
    } // switch

  std::cout << "Viewing field '" << _label << "' "<< label << ".\n"
	    << "  vector field type: " << vecFieldString << "\n"
	    << "  scale: " << _scale << "\n"
	    << "  dimensionalize flag: " << _dimensionsOkay << std::endl;
  if (!_section.isNull())
    _section->view(label);
} // view

// ----------------------------------------------------------------------
// Create PETSc vector for field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::createVector(void)
{ // createVector
  PetscErrorCode err = 0;

  if (0 != _vector) {
    err = VecDestroy(&_vector); _vector = 0;
    CHECK_PETSC_ERROR(err);
  } // if
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("GlobalOrder");

  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<typename mesh_type::SieveMesh::order_type>& order = 
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, 
					    _section->getName(), _section);
  assert(!order.isNull());

  err = VecCreate(_mesh.comm(), &_vector);
  CHECK_PETSC_ERROR(err);

  err = VecSetSizes(_vector, order->getLocalSize(), order->getGlobalSize());
  CHECK_PETSC_ERROR(err);

  err = VecSetFromOptions(_vector); CHECK_PETSC_ERROR(err);  
  logger.stagePop();
} // createVector

// ----------------------------------------------------------------------
// Create PETSc vector scatter for field. This is used to transfer
// information from the "global" PETSc vector view to the "local"
// Sieve section view.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::createScatter(void)
{ // createScatter
  assert(!_section.isNull());
  assert(!_mesh.sieveMesh().isNull());

  PetscErrorCode err = 0;
  if (0 != _scatter) {
    err = VecScatterDestroy(&_scatter); _scatter = 0;
    CHECK_PETSC_ERROR(err);
  } // if
  err = DMMeshCreateGlobalScatter(_mesh.sieveMesh(), _section, &_scatter);
  CHECK_PETSC_ERROR(err);

  if (0 != _scatterVec) {
    err = VecDestroy(&_scatterVec); _scatterVec = 0;
    CHECK_PETSC_ERROR(err);
  } // if

  // Create scatter Vec
  if (_section->sizeWithBC() > 0) {
    err = VecCreateSeqWithArray(PETSC_COMM_SELF,
				_section->getStorageSize(),
				_section->restrictSpace(),
				&_scatterVec);
    CHECK_PETSC_ERROR(err);
  } else {
    err = VecCreateSeqWithArray(PETSC_COMM_SELF, 0, 0,
				&_scatterVec);
    CHECK_PETSC_ERROR(err);
  } // else
} // createScatter

// ----------------------------------------------------------------------
// Scatter section information across processors to update the
//  PETSc vector view of the field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::scatterSectionToVector(void) const
{ // scatterSectionToVector
  scatterSectionToVector(_vector);
} // scatterSectionToVector

// ----------------------------------------------------------------------
// Scatter section information across processors to update the
//  PETSc vector view of the field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::scatterSectionToVector(const PetscVec vector) const
{ // scatterSectionToVector
  assert(!_section.isNull());
  assert(0 != _scatter);
  assert(0 != _scatterVec);
  assert(0 != vector);

  PetscErrorCode err = 0;
  err = VecScatterBegin(_scatter, _scatterVec, vector,
			INSERT_VALUES, SCATTER_FORWARD); CHECK_PETSC_ERROR(err);
  err = VecScatterEnd(_scatter, _scatterVec, vector,
		      INSERT_VALUES, SCATTER_FORWARD); CHECK_PETSC_ERROR(err);
} // scatterSectionToVector

// ----------------------------------------------------------------------
// Scatter PETSc vector information across processors to update the
// section view of the field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::scatterVectorToSection(void) const
{ // scatterVectorToSection
  scatterVectorToSection(_vector);
} // scatterVectorToSection

// ----------------------------------------------------------------------
// Scatter PETSc vector information across processors to update the
// section view of the field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::scatterVectorToSection(const PetscVec vector) const
{ // scatterVectorToSection
  assert(!_section.isNull());
  assert(0 != _scatter);
  assert(0 != _scatterVec);
  assert(0 != vector);

  PetscErrorCode err = 0;
  err = VecScatterBegin(_scatter, vector, _scatterVec,
			INSERT_VALUES, SCATTER_REVERSE); CHECK_PETSC_ERROR(err);
  err = VecScatterEnd(_scatter, vector, _scatterVec,
		      INSERT_VALUES, SCATTER_REVERSE); CHECK_PETSC_ERROR(err);
} // scatterVectorToSection

// ----------------------------------------------------------------------
// Setup split field with all one space per spatial dimension.
template<typename mesh_type>
void 
pylith::topology::Field<mesh_type>::splitDefault(void)
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
      assert(spaceDim == _section->getFiberDimension(*c_iter));
      _section->setFiberDimension(*c_iter, 1, fibration);
    } // for
} // splitDefault

// End of file 
