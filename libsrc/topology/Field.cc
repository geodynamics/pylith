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
  _vecFieldType(OTHER),
  _dimensionsOkay(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
template<typename mesh_type>
pylith::topology::Field<mesh_type>::~Field(void)
{ // destructor
  PetscErrorCode err = 0;
  if (0 != _vector) {
    err = VecDestroy(_vector); _vector = 0;
    CHECK_PETSC_ERROR(err);
  } // if

  if (0 != _scatter) {
    err = VecScatterDestroy(_scatter); _scatter = 0;
    CHECK_PETSC_ERROR(err);
  } // if
} // destructor

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
// Create seive section.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(void)
{ // newSection
  _section = new RealSection(_mesh.comm(), _mesh.debug());  
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(
				       const ALE::Obj<label_sequence>& points,
				       const int fiberDim)
{ // newSection
  typedef typename mesh_type::SieveMesh::point_type point_type;

  if (fiberDim < 0) {
    std::ostringstream msg;
    msg
      << "Fiber dimension (" << fiberDim << ") for field '" << _label
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
  } else {
    // Create empty chart
    _section->setChart(chart_type(0, 0));
    _section->setFiberDimension(points, fiberDim);  
  } // if/else
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(const DomainEnum domain,
				    const int fiberDim)
{ // newSection
  const ALE::Obj<SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  ALE::Obj<label_sequence> points;
  if (VERTICES_FIELD == domain)
    points = sieveMesh->depthStratum(0);
  else if (CELLS_FIELD == domain)
    points = sieveMesh->heightStratum(1);
  else {
    std::cerr << "Unknown value for DomainEnum: " << domain << std::endl;
    assert(0);
  } // else

  newSection(points, fiberDim);
} // newSection

// ----------------------------------------------------------------------
// Create section given chart.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(const chart_type& chart,
					       const int fiberDim)
{ // newSection
  if (_section.isNull())
    newSection();

  _section->setChart(chart);

  const typename chart_type::const_iterator chartEnd = chart.end();
  for (typename chart_type::const_iterator c_iter = chart.begin();
       c_iter != chartEnd;
       ++c_iter)
    _section->setFiberDimension(*c_iter, fiberDim);
  allocate();
} // newSection

// ----------------------------------------------------------------------
// Create section with same layout as another section.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::newSection(const Field& src)
{ // newSection
  _vecFieldType = src._vecFieldType;

  const ALE::Obj<RealSection>& srcSection = src.section();
  if (!srcSection.isNull() && _section.isNull())
    newSection();

  if (!_section.isNull()) {
    // ERROR? ARE THESE LINES CORRECT?
    _section->setAtlas(srcSection->getAtlas());
    _section->allocateStorage();
    _section->setBC(srcSection->getBC());

    if (0 != src._scatter) {
      _scatter = src._scatter;
      PetscErrorCode err = PetscObjectReference((PetscObject) _scatter);
      CHECK_PETSC_ERROR(err);
    } // if
  } // if
} // newSection

// ----------------------------------------------------------------------
// Clear variables associated with section.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::clear(void)
{ // clear
  if (!_section.isNull())
    _section->clear();

  _scale = 1.0;
  _vecFieldType = OTHER;
  _dimensionsOkay = false;
} // clear

// ----------------------------------------------------------------------
// Allocate Sieve section.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::allocate(void)
{ // allocate
  assert(!_section.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  sieveMesh->allocate(_section);
} // allocate

// ----------------------------------------------------------------------
// Zero section values.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::zero(void)
{ // zero
  if (!_section.isNull())
    _section->zero();
} // zero

// ----------------------------------------------------------------------
// Complete section by assembling across processors.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::complete(void)
{ // complete
  const ALE::Obj<SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  if (!_section.isNull())
    ALE::Completion::completeSectionAdd(sieveMesh->getSendOverlap(),
					sieveMesh->getRecvOverlap(), 
					_section, _section);
} // complete

// ----------------------------------------------------------------------
// Copy field values and metadata.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::copy(const Field& field)
{ // copy
  // Check compatibility of sections
  const int srcSize = (!field._section.isNull()) ? field._section->size() : 0;
  const int dstSize = (!_section.isNull()) ? _section->size() : 0;
  if (field.spaceDim() != spaceDim() ||
      field._vecFieldType != _vecFieldType ||
      field._scale != _scale ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot copy values from section '" << field._label 
	<< "' to section '" << _label << "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    space dim: " << field.spaceDim() << "\n"
	<< "    vector field type: " << field._vecFieldType << "\n"
	<< "    scale: " << field._scale << "\n"
	<< "    size: " << srcSize
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
    const typename chart_type::const_iterator chartEnd = chart.end();

    for (typename chart_type::const_iterator c_iter = chart.begin();
	 c_iter != chartEnd;
	 ++c_iter) {
      assert(field._section->getFiberDimension(*c_iter) ==
	     _section->getFiberDimension(*c_iter));
      _section->updatePoint(*c_iter, field._section->restrictPoint(*c_iter));
    } // for
  } // if
} // copy

// ----------------------------------------------------------------------
// Add two fields, storing the result in one of the fields.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::operator+=(const Field& field)
{ // operator+=
  // Check compatibility of sections
  const int srcSize = (!field._section.isNull()) ? field._section->size() : 0;
  const int dstSize = (!_section.isNull()) ? _section->size() : 0;
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
	<< "    size: " << srcSize
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
    const typename chart_type::const_iterator chartEnd = chart.end();

    // Assume fiber dimension is uniform
    const int fiberDim = _section->getFiberDimension(*chart.begin());
    double_array values(fiberDim);

    for (typename chart_type::const_iterator c_iter = chart.begin();
	 c_iter != chartEnd;
	 ++c_iter) {
      assert(fiberDim == field._section->getFiberDimension(*c_iter));
      assert(fiberDim == _section->getFiberDimension(*c_iter));
      field._section->restrictPoint(*c_iter, &values[0], values.size());
      _section->updateAddPoint(*c_iter, &values[0]);
    } // for
  } // if
} // operator+=

// ----------------------------------------------------------------------
// Dimensionalize field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::dimensionalize(void)
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
    const int fiberDim = _section->getFiberDimension(*chart.begin());
    double_array values(fiberDim);

    spatialdata::units::Nondimensional normalizer;

    for (typename chart_type::const_iterator c_iter = chart.begin();
	 c_iter != chartEnd;
	 ++c_iter) {
      assert(fiberDim == _section->getFiberDimension(*c_iter));
      
      _section->restrictPoint(*c_iter, &values[0], values.size());
      normalizer.dimensionalize(&values[0], values.size(), _scale);
      _section->updatePoint(*c_iter, &values[0]);
    } // for
  } // if
} // dimensionalize

// ----------------------------------------------------------------------
// Print field to standard out.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::view(const char* label)
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
    err = VecDestroy(_vector); _vector = 0;
    CHECK_PETSC_ERROR(err);
  } // if

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
    err = VecScatterDestroy(_scatter); _scatter = 0;
    CHECK_PETSC_ERROR(err);
  } // if

  err = MeshCreateGlobalScatter(_mesh.sieveMesh(), _section, &_scatter);
  CHECK_PETSC_ERROR(err);
} // createScatter

// ----------------------------------------------------------------------
// Scatter section information across processors to update the
//  PETSc vector view of the field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::scatterSectionToVector(void) const
{ // scatterSectionToVector
  assert(!_section.isNull());
  assert(0 != _scatter);
  assert(0 != _vector);

  PetscErrorCode err = 0;
  PetscVec localVec = 0;
  err = VecCreateSeqWithArray(PETSC_COMM_SELF,
			      _section->sizeWithBC(), _section->restrictSpace(),
			      &localVec); CHECK_PETSC_ERROR(err);
  err = VecScatterBegin(_scatter, localVec, _vector,
			INSERT_VALUES, SCATTER_FORWARD); CHECK_PETSC_ERROR(err);
  err = VecScatterEnd(_scatter, localVec, _vector,
		      INSERT_VALUES, SCATTER_FORWARD); CHECK_PETSC_ERROR(err);
  err = VecDestroy(localVec); CHECK_PETSC_ERROR(err);
} // scatterSectionToVector

// ----------------------------------------------------------------------
// Scatter PETSc vector information across processors to update the
// section view of the field.
template<typename mesh_type>
void
pylith::topology::Field<mesh_type>::scatterVectorToSection(void) const
{ // scatterVectorToSection
  assert(!_section.isNull());
  assert(0 != _scatter);
  assert(0 != _vector);

  PetscErrorCode err = 0;
  PetscVec localVec = 0;
  err = VecCreateSeqWithArray(PETSC_COMM_SELF,
			      _section->sizeWithBC(), _section->restrictSpace(),
			      &localVec); CHECK_PETSC_ERROR(err);
  err = VecScatterBegin(_scatter, _vector, localVec,
			INSERT_VALUES, SCATTER_REVERSE); CHECK_PETSC_ERROR(err);
  err = VecScatterEnd(_scatter, _vector, localVec,
		      INSERT_VALUES, SCATTER_REVERSE); CHECK_PETSC_ERROR(err);
  err = VecDestroy(localVec); CHECK_PETSC_ERROR(err);
} // scatterVectorToSection


// End of file 
