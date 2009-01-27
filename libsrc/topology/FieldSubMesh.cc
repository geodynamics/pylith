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

#include "FieldSubMesh.hh" // implementation of class methods

#include "SubMesh.hh" // HASA SubMesh
#include "pylith/utils/array.hh" // USES double_array

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::FieldSubMesh::FieldSubMesh(const SubMesh& mesh) :
  _mesh(mesh)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::FieldSubMesh::~FieldSubMesh(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get spatial dimension of domain.
int
pylith::topology::FieldSubMesh::spaceDim(void) const
{ // spaceDim
  const spatialdata::geocoords::CoordSys* cs = _mesh.coordsys();
  return (0 != cs) ? cs->spaceDim() : 0;
} // spaceDim

// ----------------------------------------------------------------------
// Create seive section.
void
pylith::topology::FieldSubMesh::newSection(void)
{ // newSection
  _section = new SubMeshRealSection(_mesh.comm(), _mesh.debug());  
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion.
void
pylith::topology::FieldSubMesh::newSection(
			  const ALE::Obj<SieveSubMesh::label_sequence>& points,
			  const int fiberDim)
{ // newSection
  if (fiberDim < 0) {
    std::ostringstream msg;
    msg
      << "Fiber dimension (" << fiberDim << ") for field '" << _name
      << "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  _section = new SubMeshRealSection(_mesh.comm(), _mesh.debug());

  const SieveSubMesh::point_type pointMin = 
    *std::min_element(points->begin(), points->end());
  const SieveSubMesh::point_type pointMax = 
    *std::max_element(points->begin(), points->end());
  _section->setChart(SubMeshRealSection::chart_type(pointMin, pointMax+1));
  _section->setFiberDimension(points, fiberDim);  
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion.
void
pylith::topology::FieldSubMesh::newSection(const DomainEnum domain,
					   const int fiberDim)
{ // newSection
  const ALE::Obj<SieveSubMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  ALE::Obj<SieveSubMesh::label_sequence> points;
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
void
pylith::topology::FieldSubMesh::newSection(
			const SubMeshRealSection::chart_type& chart,
			const int fiberDim)
{ // newSection
  if (_section.isNull())
    FieldSubMesh::newSection();

  _section->setChart(chart);

  const SubMeshRealSection::chart_type::const_iterator chartEnd = chart.end();
  for (SubMeshRealSection::chart_type::const_iterator c_iter = chart.begin();
       c_iter != chartEnd;
       ++c_iter)
    _section->setFiberDimension(*c_iter, fiberDim);
  allocate();
} // newSection

// ----------------------------------------------------------------------
// Create section with same layout as another section.
void
pylith::topology::FieldSubMesh::newSection(const FieldSubMesh& src)
{ // newSection
  _vecFieldType = src._vecFieldType;

  const ALE::Obj<SubMeshRealSection>& srcSection = src.section();
  if (!srcSection.isNull() && _section.isNull())
    newSection();

  if (!_section.isNull()) {
    _section->setAtlas(srcSection->getAtlas());
    _section->allocateStorage();
    _section->setBC(srcSection->getBC());
  } // if
} // newSection

// ----------------------------------------------------------------------
// Clear variables associated with section.
void
pylith::topology::FieldSubMesh::clear(void)
{ // clear
  if (!_section.isNull())
    _section->clear();

  _scale = 1.0;
  _vecFieldType = OTHER;
  _dimensionsOkay = false;
} // clear

// ----------------------------------------------------------------------
// Allocate Sieve section.
void
pylith::topology::FieldSubMesh::allocate(void)
{ // allocate
  assert(!_section.isNull());

  const ALE::Obj<SieveSubMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  sieveMesh->allocate(_section);
} // allocate

// ----------------------------------------------------------------------
// Zero section values.
void
pylith::topology::FieldSubMesh::zero(void)
{ // zero
  if (!_section.isNull())
    _section->zero();
} // zero

// ----------------------------------------------------------------------
// Complete section by assembling across processors.
void
pylith::topology::FieldSubMesh::complete(void)
{ // complete
  const ALE::Obj<SieveSubMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  if (!_section.isNull())
    ALE::Completion::completeSectionAdd(sieveMesh->getSendOverlap(),
					sieveMesh->getRecvOverlap(), 
					_section, _section);
} // complete

// ----------------------------------------------------------------------
// Copy field values and metadata.
void
pylith::topology::FieldSubMesh::copy(const FieldSubMesh& field)
{ // copy
  // Check compatibility of sections
  const int srcSize = (!field._section.isNull()) ? field._section->size() : 0;
  const int dstSize = (!_section.isNull()) ? _section->size() : 0;
  if (field.spaceDim() != spaceDim() ||
      field._vecFieldType != _vecFieldType ||
      field._scale != _scale ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot copy values from section '" << field._name 
	<< "' to section '" << _name << "'. Sections are incompatible.\n"
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
    const SubMeshRealSection::chart_type& chart = _section->getChart();
    const SubMeshRealSection::chart_type::const_iterator chartEnd = chart.end();

    for (SubMeshRealSection::chart_type::const_iterator c_iter = chart.begin();
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
void
pylith::topology::FieldSubMesh::operator+=(const FieldSubMesh& field)
{ // operator+=
  // Check compatibility of sections
  const int srcSize = (!field._section.isNull()) ? field._section->size() : 0;
  const int dstSize = (!_section.isNull()) ? _section->size() : 0;
  if (field.spaceDim() != spaceDim() ||
      field._vecFieldType != _vecFieldType ||
      field._scale != _scale ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot add values from section '" << field._name 
	<< "' to section '" << _name << "'. Sections are incompatible.\n"
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
    const SubMeshRealSection::chart_type& chart = _section->getChart();
    const SubMeshRealSection::chart_type::const_iterator chartEnd = chart.end();

    // Assume fiber dimension is uniform
    const int fiberDim = _section->getFiberDimension(*chart.begin());
    double_array values(fiberDim);

    for (SubMeshRealSection::chart_type::const_iterator c_iter = chart.begin();
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
void
pylith::topology::FieldSubMesh::dimensionalize(void)
{ // dimensionalize
  if (!_dimensionsOkay) {
    std::ostringstream msg;
    msg << "Cannot dimensionalize field '" << _name << "' because the flag "
	<< "has been set to keep field nondimensional.";
    throw std::runtime_error(msg.str());
  } // if

  if (!_section.isNull()) {
    const SubMeshRealSection::chart_type& chart = _section->getChart();
    const SubMeshRealSection::chart_type::const_iterator chartEnd = chart.end();

    // Assume fiber dimension is uniform
    const int fiberDim = _section->getFiberDimension(*chart.begin());
    double_array values(fiberDim);

    spatialdata::units::Nondimensional normalizer;

    for (SubMeshRealSection::chart_type::const_iterator c_iter = chart.begin();
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
void
pylith::topology::FieldSubMesh::view(const char* label)
{ // view
  std::string vecFieldSubMeshString;
  switch(_vecFieldType)
    { // switch
    case SCALAR:
      vecFieldSubMeshString = "scalar";
      break;
    case VECTOR:
      vecFieldSubMeshString = "vector";
      break;
    case TENSOR:
      vecFieldSubMeshString = "tensor";
      break;
    case OTHER:
      vecFieldSubMeshString = "other";
      break;
    case MULTI_SCALAR:
      vecFieldSubMeshString = "multiple scalars";
      break;
    case MULTI_VECTOR:
      vecFieldSubMeshString = "multiple vectors";
      break;
    case MULTI_TENSOR:
      vecFieldSubMeshString = "multiple tensors";
      break;
    case MULTI_OTHER:
      vecFieldSubMeshString = "multiple other values";
      break;
    default :
      std::cerr << "Unknown vector field value '" << _vecFieldType
		<< "'." << std::endl;
      assert(0);
    } // switch

  std::cout << "Viewing field '" << _name << "' "<< label << ".\n"
	    << "  vector field type: " << vecFieldSubMeshString << "\n"
	    << "  scale: " << _scale << "\n"
	    << "  dimensionalize flag: " << _dimensionsOkay << std::endl;
  if (!_section.isNull())
    _section->view(label);
} // view


// End of file 
