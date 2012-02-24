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

#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
template<typename mesh_type>
pylith::topology::FieldsNew<mesh_type>::FieldsNew(const mesh_type& mesh) :
  _mesh(mesh)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
template<typename mesh_type>
pylith::topology::FieldsNew<mesh_type>::~FieldsNew(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type>
void
pylith::topology::FieldsNew<mesh_type>::deallocate(void)
{ // deallocate
  if (!_section.isNull())
    _section->clear();

  const typename map_type::const_iterator fieldsEnd = _fields.end();
  for (typename map_type::iterator f_iter = _fields.begin(); 
       f_iter != fieldsEnd;
       ++f_iter) {
    delete f_iter->second.field; f_iter->second.field = 0;
  } // for
} // deallocate

// ----------------------------------------------------------------------
// Add field.
template<typename mesh_type>
void
pylith::topology::FieldsNew<mesh_type>::add(const char* name,
					  const char* label,
					  const int fiberDim,
					  const FieldBase::VectorFieldEnum vectorFieldType,
					  const double scale,
					  const bool dimsOkay)
{ // add
  if (hasField(name)) {
    std::ostringstream msg;
    msg << "Could not add field '" << name
	<< "' to fields manager, because it already exists.";
    throw std::runtime_error(msg.str());
  } // if

  // Set metadata
  FieldInfo info;
  info.metadata.label = label;
  info.metadata.vectorFieldType = vectorFieldType;
  info.metadata.scale = scale;
  info.metadata.dimsOkay = dimsOkay;
  
  // Set fibration and fiber dimension
  info.fibration = _fields.size();
  info.fiberDim = fiberDim;
  int sindex = 0;
  const typename map_type::const_iterator fieldsEnd = _fields.end();
  for (typename map_type::iterator f_iter=_fields.begin();
       f_iter != fieldsEnd;
       ++f_iter)
    sindex += f_iter->second.fiberDim;
  info.sindex = sindex;

  info.field = 0;

  _fields[name] = info;
} // add

// ----------------------------------------------------------------------
// Create and allocate Sieve section.
template<typename mesh_type>
void
pylith::topology::FieldsNew<mesh_type>::allocate(const ALE::Obj<typename mesh_type::SieveMesh::label_sequence>& points)
{ // allocate
  typedef typename mesh_type::SieveMesh::point_type point_type;
  typedef typename mesh_type::SieveMesh::label_sequence label_sequence;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Fields");

  // Set fiber dimension
  const int dim = fiberDim();
  assert(dim > 0);
  _section = new section_type(_mesh.comm(), _mesh.debug());
  assert(!_section.isNull());

  assert(!points.isNull());
  if (points->size() > 0) {
    const point_type pointMin = 
      *std::min_element(points->begin(), points->end());
    const point_type pointMax = 
      *std::max_element(points->begin(), points->end());
    _section->setChart(typename section_type::chart_type(pointMin, pointMax+1));
    const typename label_sequence::const_iterator pointsEnd = points->end();
    for (typename label_sequence::const_iterator p_iter = points->begin();
	 p_iter != pointsEnd;
	 ++p_iter)
      _section->setFiberDimension(*p_iter, dim);
    
  } else // Create empty chart
    _section->setChart(typename section_type::chart_type(0, 0));

  // Set spaces
  const typename map_type::const_iterator fieldsEnd = _fields.end();
  for (typename map_type::iterator f_iter=_fields.begin();
       f_iter != fieldsEnd;
       ++f_iter)
    _section->addSpace();
  for (typename map_type::iterator f_iter=_fields.begin();
       f_iter != fieldsEnd;
       ++f_iter)
    _section->spaceFiberDimension(f_iter->second.fibration, 
				  f_iter->second.fiberDim);
  // Allocate section
  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  sieveMesh->allocate(_section);

  logger.stagePop();
} // allocate

// ----------------------------------------------------------------------
// Create and allocate Sieve section.
template<typename mesh_type>
void
pylith::topology::FieldsNew<mesh_type>::allocate(const int_array& points)
{ // allocate
  typedef typename mesh_type::SieveMesh::point_type point_type;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Field");

  // Set fiber dimension
  const int dim = fiberDim();
  assert(dim > 0);
  _section = new section_type(_mesh.comm(), _mesh.debug());
  assert(!_section.isNull());

  const int npts = points.size();
  if (npts > 0) {
    const point_type pointMin = points.min();
    const point_type pointMax = points.max();
    _section->setChart(typename section_type::chart_type(pointMin, pointMax+1));
    for (int i=0; i < npts; ++i)
      _section->setFiberDimension(points[i], dim);

  } else  // create empty chart
    _section->setChart(typename section_type::chart_type(0, 0));

  // Set spaces
  const typename map_type::const_iterator fieldsEnd = _fields.end();
  for (typename map_type::iterator f_iter=_fields.begin();
       f_iter != fieldsEnd;
       ++f_iter)
    _section->addSpace();
  for (typename map_type::iterator f_iter=_fields.begin();
       f_iter != fieldsEnd;
       ++f_iter)
    _section->spaceFiberDimension(f_iter->second.fibration, 
				  f_iter->second.fiberDim);

  // Allocate section
  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  sieveMesh->allocate(_section);

  logger.stagePop();
} // allocate

// ----------------------------------------------------------------------
// Create and allocate Sieve section.
template<typename mesh_type>
void
pylith::topology::FieldsNew<mesh_type>::allocate(const FieldBase::DomainEnum domain,
					       const int stratum)
{ // allocate
  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  ALE::Obj<typename mesh_type::SieveMesh::label_sequence> points;
  if (FieldBase::VERTICES_FIELD == domain)
    points = sieveMesh->depthStratum(stratum);
  else if (FieldBase::CELLS_FIELD == domain)
    points = sieveMesh->heightStratum(stratum);
  else {
    std::cerr << "Unknown value for DomainEnum: " << domain << std::endl;
    assert(0);
    throw std::logic_error("Bad domain enum in Field.");
  } // else

  allocate(points);
} // allocate

// ----------------------------------------------------------------------
// Get field.
template<typename mesh_type>
pylith::topology::Field<mesh_type>&
pylith::topology::FieldsNew<mesh_type>::get(const char* name)
{ // get
  typename map_type::iterator f_iter = _fields.find(name);
  if (f_iter == _fields.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name
	<< "' in fields manager for retrieval.";
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
// Compute total fiber dimension for section.
template<typename mesh_type>
int
pylith::topology::FieldsNew<mesh_type>::fiberDim(void) const
{ // fiberDim
  int fiberDim = 0;
  const typename map_type::const_iterator fieldsEnd = _fields.end();
  for (typename map_type::const_iterator f_iter=_fields.begin();
       f_iter != fieldsEnd;
       ++f_iter)
    fiberDim += f_iter->second.fiberDim;

  if (fiberDim < 0) {
    std::ostringstream msg;
    msg << "Fiber dimension (" << fiberDim << ") for Fields object must "
	<< "be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if


  return fiberDim;
} // fiberDim

// ----------------------------------------------------------------------
// Get index of first value of field in section.
template<typename mesh_type>
int
pylith::topology::FieldsNew<mesh_type>::sectionIndex(const char* name) const
{ // sectionIndex
  typename map_type::const_iterator f_iter = _fields.find(name);
  if (f_iter == _fields.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name
	<< "' in fields manager for retrieval of section index.";
    throw std::runtime_error(msg.str());
  } // if

  return f_iter->second.sindex;
} // sectionIndex

// ----------------------------------------------------------------------
// Get fiber dimension of field in section.
template<typename mesh_type>
int
pylith::topology::FieldsNew<mesh_type>::sectionFiberDim(const char* name) const
{ // sectionFiberDim
  typename map_type::const_iterator f_iter = _fields.find(name);
  if (f_iter == _fields.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name
	<< "' in fields manager for retrieval of section index.";
    throw std::runtime_error(msg.str());
  } // if

  return f_iter->second.fiberDim;
} // sectionFiberDim

// ----------------------------------------------------------------------
// Complete section by assembling across processors.
template<typename mesh_type>
void
pylith::topology::FieldsNew<mesh_type>::complete(void)
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
// Get names of all fields
template<typename mesh_type>
void
pylith::topology::FieldsNew<mesh_type>::fieldNames(int *numNames, 
						   char*** names) const
{ // fieldNames
  assert(numNames);
  assert(names);

  *numNames = _fields.size();
  *names = new char*[_fields.size()];
  assert(*names);
  const typename map_type::const_iterator namesEnd = _fields.end();
  int i = 0;
  for (typename map_type::const_iterator name = _fields.begin(); 
       name != namesEnd;
       ++name) {
    const char len = name->first.length();
    char* newName = 0;
    if (len > 0) {
      newName = new char[len+1];
      strncpy(newName, name->first.c_str(), len+1);
    } else {
      newName = new char[1];
      newName[0] ='\0';
    } // if/else
    (*names)[i++] = newName;
  } // for
} // fieldNames

// ----------------------------------------------------------------------
// View fields and section.
template<typename mesh_type>
void
pylith::topology::FieldsNew<mesh_type>::view(const char* label)
{ // view
  std::cout << "Fields '" << label << "':\n";
  const typename map_type::const_iterator fieldsEnd = _fields.end();
  for (typename map_type::const_iterator f_iter = _fields.begin();
       f_iter != fieldsEnd;
       ++f_iter)
    std::cout << "  Field: " << f_iter->first
	      << ", fibration: " << f_iter->second.fibration
	      << ", fiber dim: " << f_iter->second.fiberDim
	      << ", first value index: " << f_iter->second.sindex
	      << std::endl;
  _section->view("Section");
} // view


// End of file 
