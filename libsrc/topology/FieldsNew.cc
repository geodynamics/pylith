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

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Fields");

  // Set fiber dimension
  const int dim = fiberDim();
  assert(dim > 0);
#if defined(USE_UNIFORMSECTION)
  _section = new section_type(_mesh.comm(), dim, _mesh.debug());
#else
  _section = new section_type(_mesh.comm(), _mesh.debug());
#endif
  assert(!_section.isNull());

  // Set spaces
  const typename map_type::const_iterator fieldsEnd = _fields.end();
  for (typename map_type::iterator f_iter=_fields.begin();
       f_iter != fieldsEnd;
       ++f_iter)
    _section->addSpace();

  assert(!points.isNull());
  if (points->size() > 0) {
    const point_type pointMin = 
      *std::min_element(points->begin(), points->end());
    const point_type pointMax = 
      *std::max_element(points->begin(), points->end());
    _section->setChart(typename section_type::chart_type(pointMin, pointMax+1));
    _section->setFiberDimension(points, dim);
    
    int fibration = 0;
    for (typename map_type::iterator f_iter=_fields.begin();
	 f_iter != fieldsEnd;
	 ++f_iter, ++fibration)
      _section->setFiberDimension(points, f_iter->second.fiberDim, fibration);
  } else // Create empty chart
    _section->setChart(typename section_type::chart_type(0, 0));

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
#if defined(USE_UNIFORMSECTION)
  _section = new section_type(_mesh.comm(), dim, _mesh.debug());
#else
  _section = new section_type(_mesh.comm(), _mesh.debug());
#endif
  assert(!_section.isNull());

  // Set spaces
  const typename map_type::const_iterator fieldsEnd = _fields.end();
  for (typename map_type::iterator f_iter=_fields.begin();
       f_iter != fieldsEnd;
       ++f_iter)
    _section->addSpace();

  const int npts = points.size();
  if (npts > 0) {
    const point_type pointMin = points.min();
    const point_type pointMax = points.max();
    _section->setChart(typename section_type::chart_type(pointMin, pointMax+1));
    for (int i=0; i < npts; ++i)
      _section->setFiberDimension(points[i], dim);

    int fibration = 0;
    for (typename map_type::iterator f_iter=_fields.begin();
	 f_iter != fieldsEnd;
	 ++f_iter, ++fibration)
      for (int i=0; i < npts; ++i)
	_section->setFiberDimension(points[i], 
				    f_iter->second.fiberDim, fibration);
  } else  // create empty chart
    _section->setChart(typename section_type::chart_type(0, 0));

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
const pylith::topology::Field<mesh_type>&
pylith::topology::FieldsNew<mesh_type>::get(const char* name) const
{ // get
  typename map_type::const_iterator f_iter = _fields.find(name);
  if (f_iter == _fields.end()) {
    std::ostringstream msg;
    msg << "Could not find field '" << name
	<< "' in fields manager for retrieval.";
    throw std::runtime_error(msg.str());
  } // if

  const int fibration = f_iter->second.fibration;
  assert(fibration >= 0 && fibration < _fields.size());

  delete f_iter->second.field; f_iter->second.field = 0;
  assert(!_section.isNull());
  f_iter->second.field = 
    new Field<mesh_type>(_mesh, _section->getFibration(fibration), 
			 f_iter->second.metadata);
  assert(0 != f_iter->second.field);

  return *f_iter->second.field;
} // get
	   
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

  delete f_iter->second.field; f_iter->second.field = 0;
  assert(!_section.isNull());
  f_iter->second.field = 
    new Field<mesh_type>(_mesh, _section->getFibration(fibration), 
			 f_iter->second.metadata);
  assert(0 != f_iter->second.field);

  return *f_iter->second.field;
} // get

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
// Get names of all fields
template<typename mesh_type>
void
pylith::topology::FieldsNew<mesh_type>::fieldNames(int *numNames, 
						   std::string** names) const
{ // fieldNames
  const int size = _fields.size();
  *numNames = size;
  *names = new std::string[size];
  
  int i = 0;
  const typename map_type::const_iterator fieldsEnd = _fields.end();
  for (typename map_type::const_iterator f_iter = _fields.begin();
       f_iter != fieldsEnd;
       ++f_iter)
    (*names)[i++] = f_iter->first;
} // fieldNames

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


// End of file 
