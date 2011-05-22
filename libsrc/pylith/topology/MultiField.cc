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

#include "MultiField.hh" // implementation of class methods

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
pylith::topology::MultiField<mesh_type, section_type>::MultiField(const mesh_type& mesh) :
  _mesh(mesh)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
template<typename mesh_type, typename section_type>
pylith::topology::MultiField<mesh_type, section_type>::~MultiField(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiField<mesh_type, section_type>::deallocate(void)
{ // deallocate
  PetscErrorCode err = 0;
  
  if (!_section.isNull())
    _section->clear();

  const typename scatter_map_type::const_iterator scattersEnd = _scatters.end();
  for (typename scatter_map_type::iterator s_iter=_scatters.begin();
       s_iter != scattersEnd;
       ++s_iter) {
    if (s_iter->second.vector) {
      err = VecDestroy(&s_iter->second.vector);CHECK_PETSC_ERROR(err);
    } // if

    if (s_iter->second.scatter) {
      err = VecScatterDestroy(&s_iter->second.scatter);CHECK_PETSC_ERROR(err);
    } // if

    if (s_iter->second.scatterVec) {
      err = VecDestroy(&s_iter->second.scatterVec);CHECK_PETSC_ERROR(err);
    } // if
  } // for
  _scatters.clear();
} // deallocate

// ----------------------------------------------------------------------
// Set label of section.
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiField<mesh_type, section_type>::label(const char* value)
{ // label
  assert(!_section.isNull());
  _section->setName(value);
} // label

// ----------------------------------------------------------------------
// Get the chart size.
template<typename mesh_type, typename section_type>
int
pylith::topology::MultiField<mesh_type, section_type>::chartSize(void) const
{ // chartSize
  return _section.isNull() ? 0 : _section->getChart().size();
} // chartSize

// ----------------------------------------------------------------------
// Get the number of degrees of freedom.
template<typename mesh_type, typename section_type>
int
pylith::topology::MultiField<mesh_type, section_type>::sectionSize(void) const
{ // sectionSize
  return _section.isNull() ? 0 : _section->size();
} // sectionSize

// ----------------------------------------------------------------------
// Add field.
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiFields<mesh_type>::add(const char* name,
					      const char* label,
					      const FieldBase::VectorFieldEnum vectorFieldType,
					      const double scale,
					      const bool dimsOkay)
{ // add
  if (hasField(name)) {
    std::ostringstream msg;
    msg << "Could not add field '" << name
	<< "' to multiple fields object, because it already exists.";
    throw std::runtime_error(msg.str());
  } // if

  // Set metadata
  FieldInfo info;
  info.metadata.label = label;
  info.metadata.vectorFieldType = vectorFieldType;
  info.metadata.scale = scale;
  info.metadata.dimsOkay = dimsOkay;
  
  // Set field index.
  info.fieldIndex = _fields.size();
  info.field = 0;

  _fields[name] = info;
} // add

// ----------------------------------------------------------------------
// Get field.
template<typename mesh_type, typename section_type>
pylith::topology::Field<mesh_type>&
pylith::topology::MultiFields<mesh_type>::get(const char* name)
{ // get
  typename map_type::iterator f_iter = _fields.find(name);
  if (f_iter == _fields.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name
	<< "' in multiple fields object for retrieval.";
    throw std::runtime_error(msg.str());
  } // if
  const int fibration = f_iter->second.fibration;
  assert(fibration >= 0 && fibration < _fields.size());

  if (!f_iter->second.field) {
    delete f_iter->second.field; f_iter->second.field = 0;
    assert(!_section.isNull());
    f_iter->second.field = 
      new Field<mesh_type>(_mesh, _section->getFibration(fibration), 
			   f_iter->second.metadata);
    assert(0 != f_iter->second.field);
  } // if

  return *f_iter->second.field;
} // get

// ----------------------------------------------------------------------
// Get index of field in collection of fields.
template<typename mesh_type, typename section_type>
int
pylith::topology::MultiFields<mesh_type>::fieldIndex(const char* name) const
{ // fieldIndex
  const typename map_type::const_iterator iterTarget = _fields.find(name);
  if (iterTarget == _fields.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name
	<< "' in multiple fields object for retrieval of field index.";
    throw std::runtime_error(msg.str());
  } // if

  assert(f_iter);
  const int index = f_iter->second.fibration;
  assert(index >= 0);

  return sindex;
} // fieldIndex

// ----------------------------------------------------------------------
// Get index of first value of field in field.
template<typename mesh_type, typename section_type>
int
pylith::topology::MultiFields<mesh_type>::fieldStartIndex(const int fieldIndex,
							  const point_type point) const
{ // fieldStartIndex
  assert(_section.numSpaces() > 0);
  assert(fieldIndex >= 0 && fieldIndex < _section.numSpaces());
  
  int sindex = 0;
  for (int i=0; i < fieldIndex; ++i) {
    const ALE::Obj<RealSection>& subsection =
      _section->getFibration(i);
    sindex += subsection->getFiberDimension(point);
  } // for

  return sindex;
} // fieldStartIndex

// ----------------------------------------------------------------------
// Get fiber dimension of field in section.
template<typename mesh_type, typename section_type>
int
pylith::topology::MultiFields<mesh_type>::fieldFiberDim(const int fieldIndex,
							const point_type point) const
{ // fieldFiberDim
  assert(_section.numSpaces() > 0);
  assert(fieldIndex >= 0 && fieldIndex < _section.numSpaces());

  const ALE::Obj<RealSection>& subsection =
    _section->getFibration(fieldIndex);
  const int fiberDim = subsection->getFiberDimension(point);

  return fiberDim;
} // fieldFiberDim

// ----------------------------------------------------------------------
// Compute total fiber dimension for section.
template<typename mesh_type, typename section_type>
int
pylith::topology::MultiFields<mesh_type>::fiberDim(void) const
{ // fiberDim
  const int numSpaces = _section.numSpaces();
  assert(numSpaces > 0);

  int fiberDim = 0;
  for (int i=0; i < numSpaces; ++i) {
    const ALE::Obj<RealSection>& subsection =
      _section->getFibration(i);
    fiberDim += subsection->getFiberDimension(point);
  } // for

  return fiberDim;
} // fiberDim

// ----------------------------------------------------------------------
// Clear variables associated with section.
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiField<mesh_type, section_type>::clear(void)
{ // clear
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");

  deallocate();
  if (!_section.isNull())
    _section->clear();
  _fields.clear();

  logger.stagePop();
} // clear

// ----------------------------------------------------------------------
// Allocate Sieve section.
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiField<mesh_type, section_type>::allocate(void)
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
// Create section with same layout as another section.
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiField<mesh_type, section_type>::cloneSection(const MultiField& src)
{ // cloneSection
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");

  deallocate();

  // Copy metadata
  _fields = src._fields;

  // Section
  const ALE::Obj<section_type>& srcSection = src.section();
  if (!srcSection.isNull() && _section.isNull()) {
    logger.stagePop();
    newSection();
    logger.stagePush("Field");
  } // if

  if (!_section.isNull()) {
    // Note: We retain the original label of the section.

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
    const char* sectionLabel = _section->getName().c_str();
    if (src._scatters.size() > 0) {
      const typename scatter_map_type::const_iterator scattersEnd = src._scatters.end();
      for (typename scatter_map_type::const_iterator s_iter=src._scatters.begin();
	   s_iter != scattersEnd;
	   ++s_iter) {
	ScatterInfo sinfo;
	sinfo.vector = 0;
	sinfo.scatterVec = 0;

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
	
	err = PetscObjectSetName((PetscObject)sinfo.vector, sectionLabel);
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
// Complete section by assembling across processors.
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiField<mesh_type, section_type>::complete(void)
{ // complete
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Completion");

  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  if (!_section.isNull())
    ALE::Completion::completeSectionAdd(sieveMesh->getSendOverlap(),
					sieveMesh->getRecvOverlap(), 
					_section, _section);

  logger.stagePop();
} // complete

// ----------------------------------------------------------------------
// Zero section values (excluding constrained DOF).
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiField<mesh_type, section_type>::zero(void)
{ // zero
  if (!_section.isNull())
    _section->zero(); // Does not zero BC.
} // zero

// ----------------------------------------------------------------------
// Zero section values (including constrained DOF).
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiField<mesh_type, section_type>::zeroAll(void)
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
// Copy field values and metadata.
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiField<mesh_type, section_type>::copy(const MultiField& field)
{ // copy
  // Check compatibility of sections
  const int srcChartSize = field.chartSize();
  const int dstChartSize = chartSize();
  const int srcSectionSize = field.sectionSize();
  const int dstSectionSize = sectionSize();
  if (srcChartSize != dstChartSize ||
      srcSectionSize != dstSectionSize) {
    std::ostringstream msg;

    if (field._section.isNull())
      msg << "Cannot copy values from null field to ";
    else
      msg << "Cannot copy values from field '" << field._section->getName()
	  << "' to ";
    if (_section.isNull())
      msg << "null section.";
    else
      msg << "section '" << _section->getName() << "'.";
    msg << " Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    chart size: " << srcChartSize << "\n"
	<< "    section size: " << srcSectionSize << "\n"
	<< "  Destination section:\n"
	<< "    chart size: " << dstChartSize << "\n"
	<< "    section size: " << dstSectionSize << "\n";
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
pylith::topology::MultiField<mesh_type, section_type>::copy(const ALE::Obj<section_type>& osection)
{ // copy
  // Check compatibility of sections
  const int srcChartSize = 
    (!osection.isNull()) ? osection->getChart().size() : 0;
  const int srcSectionSize = 
    (!osection.isNull()) ? osection->size() : 0;

  const int dstChartSize = chartSize();
  const int dstSectionSize = sectionSize();
  if (srcChartSize != dstChartSize ||
      srcSectionSize != dstSectionSize) {
    std::ostringstream msg;

    if (field._section.isNull())
      msg << "Cannot copy values from null section to ";
    else
      msg << "Cannot copy values from section '" << field._section->getName()
	  << "' to ";
    if (_section.isNull())
      msg << "null section.";
    else
      msg << "section '" << _section->getName() << "'.";
    msg << " Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    chart size: " << srcChartSize << "\n"
	<< "    section size: " << srcSectionSize << "\n"
	<< "  Destination section:\n"
	<< "    chart size: " << dstChartSize << "\n"
	<< "    section size: " << dstSectionSize << "\n";
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
pylith::topology::MultiField<mesh_type, section_type>&
pylith::topology::MultiField<mesh_type, section_type>::operator+=(const MultiField& field)
{ // operator+=
  // Check compatibility of sections
  const int srcChartSize = field.chartSize();
  const int dstChartSize = chartSize();
  const int srcSectionSize = field.sectionSize();
  const int dstSectionSize = sectionSize();
  if (srcChartSize != dstChartSize ||
      srcSectionSize != dstSectionSize) {
    std::ostringstream msg;

    if (field._section.isNull())
      msg << "Cannot add values from null field to ";
    else
      msg << "Cannot add values from field '" << field._section->getName()
	  << "' to ";
    if (_section.isNull())
      msg << "null section.";
    else
      msg << "section '" << _section->getName() << "'.";
    msg << " Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    chart size: " << srcChartSize << "\n"
	<< "    section size: " << srcSectionSize << "\n"
	<< "  Destination section:\n"
	<< "    chart size: " << dstChartSize << "\n"
	<< "    section size: " << dstSectionSize << "\n";
    throw std::runtime_error(msg.str());
  } // if
  assert( (_section.isNull() && field._section.isNull()) ||
	  (!_section.isNull() && !field._section.isNull()) );

  if (!_section.isNull()) {
    // Add values from field
    const chart_type& chart = _section->getChart();
    const typename chart_type::const_iterator chartBegin = chart.begin();
    const typename chart_type::const_iterator chartEnd = chart.end();

    for (typename chart_type::const_iterator c_iter = chartBegin;
	 c_iter != chartEnd;
	 ++c_iter) {
      assert(osection->getFiberDimension(*c_iter) ==
	     _section->getFiberDimension(*c_iter));
      if (field._section->getFiberDimension(*c_iter) > 0) {
	assert(fiberDim == field._section->getFiberDimension(*c_iter));
	assert(fiberDim == _section->getFiberDimension(*c_iter));
	_section->updatePointAllAdd(*c_iter, 
				    field._section->restrictPoint(*c_iter));
      } // if
    } // for
  } // if

  return *this;
} // operator+=

// ----------------------------------------------------------------------
// Dimensionalize field.
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiField<mesh_type, section_type>::dimensionalize(void) const
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
    // :FIX THIS:
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
pylith::topology::MultiField<mesh_type, section_type>::view(const char* label) const
{ // view
  if (_section.isNull())
    std::cout << "Fields in collection of fields '" << _section->getName()
	      << "':\n";
  else
    std::cout << "Fields in unknown collection of fields:\n";

  std::string vecFieldString = "";
  const typename map_type::const_iterator fieldsEnd = _fields.end();
  for (typename map_type::const_iterator f_iter = _fields.begin();
       f_iter != fieldsEnd;
       ++f_iter) {
    std::cout << "  Field: " << f_iter->first
	      << ", index: " << f_iter->second.fieldIndex;
    
    const FieldBase::Metadata& metadata = f_iter->second.metadata;
    switch(metadata.vectorFieldType)
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
	std::cerr << "Unknown vector field value '" 
		  << metadata.vectorFieldType
		  << "'." << std::endl;
	assert(0);
	throw std::logic_error("Bad vector field type in MultiField.");
      } // switch

    std::cout << ",  vector field type: " << vecFieldString
	      << ", scale: " << metadata.scale
	      << ", dimensionalize flag: " << metadata.dimsOkay
	      << std::endl;
  _section->view("Section");
} // view

// ----------------------------------------------------------------------
// Create PETSc vector scatter for field. This is used to transfer
// information from the "global" PETSc vector view to the "local"
// Sieve section view.
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiField<mesh_type, section_type>::createScatter(const char* context)
{ // createScatter
  assert(context);
  assert(!_section.isNull());
  assert(!_mesh.sieveMesh().isNull());

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

  // Get global order (create if necessary).
  const std::string& orderLabel = _section->getName();
  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh =
    _mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<typename mesh_type::SieveMesh::order_type>& order = 
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, orderLabel,
					    _section);
  assert(!order.isNull());

  // Create scatter
  err = DMMeshCreateGlobalScatter(_mesh.sieveMesh(), _section, order, false, &sinfo.scatter);
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
void
pylith::topology::MultiField<mesh_type, section_type>::createScatter(const typename ALE::Obj<typename SieveMesh::numbering_type> numbering,
								const char* context)
{ // createScatter
  assert(!numbering.isNull());
  assert(context);
  assert(!_section.isNull());
  assert(!_mesh.sieveMesh().isNull());

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

  // Get global order (create if necessary).
  const std::string& orderLabel = 
    (strlen(context) > 0) ?
    _section->getName() + std::string("_") + std::string(context) :
    _section->getName();
  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh =
    _mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<typename mesh_type::SieveMesh::order_type>& order = 
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, orderLabel,
                                            numbering->getChart().begin(),
                                            numbering->getChart().end(),
                                            _section);
  assert(!order.isNull());

  // Create scatter
  err = DMMeshCreateGlobalScatter(_mesh.sieveMesh(), _section, order, false, &sinfo.scatter); 
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
void
pylith::topology::MultiField<mesh_type, section_type>::createScatterWithBC(const char* context)
{ // createScatterWithBC
  assert(context);
  assert(!_section.isNull());
  assert(!_mesh.sieveMesh().isNull());


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

  // Get global order (create if necessary).
  const std::string& orderLabel = _section->getName();
  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh =
    _mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<typename mesh_type::SieveMesh::order_type>& order = 
    sieveMesh->getFactory()->getGlobalOrderWithBC(sieveMesh, orderLabel,
						  _section);
  assert(!order.isNull());

  // Create scatter
  err = DMMeshCreateGlobalScatter(_mesh.sieveMesh(), _section, order, true, &sinfo.scatter); 
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
} // createScatterWithBC

// ----------------------------------------------------------------------
// Create PETSc vector scatter for field. This is used to transfer
// information from the "global" PETSc vector view to the "local"
// Sieve section view. The PETSc vector includes constrained DOF. Use
// createScatter() if constrained DOF should be omitted from the PETSc
// vector.
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiField<mesh_type, section_type>::createScatterWithBC(const typename ALE::Obj<typename SieveMesh::numbering_type> numbering,
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

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("GlobalOrder");

  // Get global order (create if necessary).
  const std::string& orderLabel = 
    (strlen(context) > 0) ?
    _section->getName() + std::string("_") + std::string(context) :
    _section->getName();
  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh =
    _mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<typename mesh_type::SieveMesh::order_type>& order = 
    sieveMesh->getFactory()->getGlobalOrderWithBC(sieveMesh, orderLabel,
                                                  numbering->getChart().begin(),
                                                  numbering->getChart().end(),
                                                  _section);
  assert(!order.isNull());
  //order->view("GLOBAL ORDER"); // DEBUG

  // Create scatter
  err = DMMeshCreateGlobalScatter(_mesh.sieveMesh(), _section, order, true, &sinfo.scatter); 
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

#if 0
  std::cout << "CONTEXT: " << context 
	    << ", orderLabel: " << orderLabel
	    << ", section size w/BC: " << _section->sizeWithBC()
	    << ", section size: " << _section->size()
	    << ", section storage size: " << _section->getStorageSize()
	    << ", global numbering size: " << numbering->getGlobalSize()
	    << ", global size: " << order->getGlobalSize()
	    << ", scatter from size: " << sinfo.scatter->from_n
	    << std::endl;
#endif
  
  logger.stagePop();
} // createScatterWithBC

// ----------------------------------------------------------------------
// Get PETSc vector associated with field.
template<typename mesh_type, typename section_type>
PetscVec
pylith::topology::MultiField<mesh_type, section_type>::vector(const char* context)
{ // vector
  ScatterInfo& sinfo = _getScatter(context);
  return sinfo.vector;
} // vector

// ----------------------------------------------------------------------
// Get PETSc vector associated with field.
template<typename mesh_type, typename section_type>
const PetscVec
pylith::topology::MultiField<mesh_type, section_type>::vector(const char* context) const
{ // vector
  const ScatterInfo& sinfo = _getScatter(context);
  return sinfo.vector;
} // vector

// ----------------------------------------------------------------------
// Scatter section information across processors to update the
//  PETSc vector view of the field.
template<typename mesh_type, typename section_type>
void
pylith::topology::MultiField<mesh_type, section_type>::scatterSectionToVector(const char* context) const
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
pylith::topology::MultiField<mesh_type, section_type>::scatterSectionToVector(const PetscVec vector,
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
pylith::topology::MultiField<mesh_type, section_type>::scatterVectorToSection(const char* context) const
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
pylith::topology::MultiField<mesh_type, section_type>::scatterVectorToSection(const PetscVec vector,
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
pylith::topology::MultiField<mesh_type, section_type>::splitDefault(void)
{ // splitDefault
#if 0
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
#else
  // :TODO: tell solver to split vector fields into pieces
  throw std::logic_error("MultiField::splitDefault() needs updating.");
#endif
} // splitDefault

// ----------------------------------------------------------------------
// Get scatter for given context.
template<typename mesh_type, typename section_type>
typename pylith::topology::MultiField<mesh_type, section_type>::ScatterInfo&
pylith::topology::MultiField<mesh_type, section_type>::_getScatter(const char* context,
							      const bool createOk)
{ // _getScatter
  assert(context);

  const bool isNewScatter = _scatters.find(context) == _scatters.end();

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
const typename pylith::topology::MultiField<mesh_type, section_type>::ScatterInfo&
pylith::topology::MultiField<mesh_type, section_type>::_getScatter(const char* context) const
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
