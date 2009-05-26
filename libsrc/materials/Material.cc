// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "Material.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/utils/array.hh" // USES double_array, std::vector

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES memcpy()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::RestrictVisitor RestrictVisitor;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::Material::Material(const int dimension,
				      const int tensorSize,
				      const Metadata& metadata) :
  _dt(0.0),
  _properties(0),
  _stateVars(0),
  _normalizer(new spatialdata::units::Nondimensional),
  _numPropsQuadPt(0),
  _numVarsQuadPt(0),
  _dimension(dimension),
  _tensorSize(tensorSize),
  _needNewJacobian(false),
  _dbProperties(0),
  _dbInitialState(0),
  _id(0),
  _label(""),
  _metadata(metadata)
{ // constructor
  const string_vector& properties = metadata.properties();
  const int numProperties = properties.size();
  for (int i=0; i < numProperties; ++i)
    _numPropsQuadPt += metadata.fiberDim(properties[i].c_str(),
					 Metadata::PROPERTY);
  assert(_numPropsQuadPt >= 0);

  const string_vector& stateVars = metadata.stateVars();
  const int numStateVars = stateVars.size();
  for (int i=0; i < numStateVars; ++i)
    _numVarsQuadPt += metadata.fiberDim(stateVars[i].c_str(),
					Metadata::STATEVAR);
  assert(_numVarsQuadPt >= 0);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::Material::~Material(void)
{ // destructor
  delete _normalizer; _normalizer = 0;
  delete _properties; _properties = 0;
  delete _stateVars; _stateVars = 0;

  // Python db object owns databases, so just set pointer to null
  // :KLUDGE: Should use shared pointer
  _dbProperties = 0;
  _dbInitialState = 0;
} // destructor

// ----------------------------------------------------------------------
// Set scales used to nondimensionalize physical properties.
void
pylith::materials::Material::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
  if (0 == _normalizer)
    _normalizer = new spatialdata::units::Nondimensional(dim);
  else
    *_normalizer = dim;
} // normalizer

// ----------------------------------------------------------------------
// Get physical property parameters and initial state (if used) from database.
void
pylith::materials::Material::initialize(
		     const topology::Mesh& mesh,
		     feassemble::Quadrature<topology::Mesh>* quadrature)
{ // initialize
  assert(0 != _dbProperties);
  assert(0 != quadrature);
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Materials");

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int numBasis = quadrature->numBasis();
  const int spaceDim = quadrature->spaceDim();

  // Get cells associated with material
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", _id);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);

  // Create field to hold physical properties.
  delete _properties; _properties = new topology::Field<topology::Mesh>(mesh);
  _properties->label("properties");
  assert(0 != _properties);
  int fiberDim = numQuadPts * _numPropsQuadPt;
  _properties->newSection(cells, fiberDim);
  _properties->allocate();
  _properties->zero();
  std::cout << "Initializing material " << _label << ": Created properties" << std::endl;
  const ALE::Obj<RealSection>& propertiesSection = _properties->section();
  assert(!propertiesSection.isNull());

#if !defined(PRECOMPUTE_GEOMETRY)
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  // Create arrays for querying.
  const int numDBProperties = _metadata.numDBProperties();
  double_array quadPtsGlobal(numQuadPts*spaceDim);
  double_array propertiesQuery(numDBProperties);
  double_array propertiesCell(numQuadPts*numDBProperties);

  // Setup database for quering for physical properties
  assert(0 != _dbProperties);
  _dbProperties->open();
  _dbProperties->queryVals(_metadata.dbProperties(),
			   _metadata.numDBProperties());

  // Create field to hold state variables. We create the field even
  // if there is no initial state, because this we will use this field
  // to hold the state variables.
  delete _stateVars; _stateVars = new topology::Field<topology::Mesh>(mesh);
  _stateVars->label("state variables");
  fiberDim = numQuadPts * _numVarsQuadPt;
  if (fiberDim > 0) {
    assert(0 != _stateVars);
    const ALE::Obj<RealSection::chart_type>& chart = 
      propertiesSection->getChart();
    assert(!chart.isNull());
    _stateVars->newSection(*chart, fiberDim);
    _stateVars->allocate();
    _stateVars->zero();
  } // if
  std::cout << "Initializing material " << _label << ": Created state vars" << std::endl;
  const ALE::Obj<RealSection>& stateVarsSection = 
    (fiberDim > 0) ? _stateVars->section() : 0;

  // Create arrays for querying
  const int numDBStateVars = _metadata.numDBStateVars();
  double_array stateVarsQuery;
  double_array stateVarsCell;
  if (0 != _dbInitialState) {
    assert(!stateVarsSection.isNull());
    assert(numDBStateVars > 0);
    assert(_numVarsQuadPt > 0);
    stateVarsQuery.resize(numDBStateVars);
    stateVarsCell.resize(numQuadPts*numDBStateVars);
    
    // Setup database for querying for initial state variables
    _dbInitialState->open();
    _dbInitialState->queryVals(_metadata.dbStateVars(),
			       _metadata.numDBStateVars());
  } // if

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
    
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
    quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    const double_array& quadPtsNonDim = quadrature->quadPts();
    quadPtsGlobal = quadPtsNonDim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				lengthScale);

    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, index=0; 
	 iQuadPt < numQuadPts; 
	 ++iQuadPt, index+=spaceDim) {
      int err = _dbProperties->query(&propertiesQuery[0], numDBProperties,
				     &quadPtsGlobal[index], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find parameters for physical properties at \n"
	    << "(";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPtsGlobal[index+i];
	msg << ") in material " << _label << "\n"
	    << "using spatial database '" << _dbProperties->label() << "'.";
	throw std::runtime_error(msg.str());
      } // if
      _dbToProperties(&propertiesCell[iQuadPt*_numPropsQuadPt], 
		      propertiesQuery);
      _nondimProperties(&propertiesCell[iQuadPt*_numPropsQuadPt],
			_numPropsQuadPt);

      if (0 != _dbInitialState) {
	err = _dbInitialState->query(&stateVarsQuery[0], numDBStateVars,
				     &quadPtsGlobal[index], spaceDim, cs);
	if (err) {
	  std::ostringstream msg;
	  msg << "Could not find initial state variables at \n" << "(";
	  for (int i=0; i < spaceDim; ++i)
	    msg << "  " << quadPtsGlobal[index+i];
	  msg << ") in material " << _label << "\n"
	      << "using spatial database '" << _dbInitialState->label()
	      << "'.";
	  throw std::runtime_error(msg.str());
	} // if
	_dbToStateVars(&stateVarsCell[iQuadPt*_numVarsQuadPt], 
		       stateVarsQuery);
	_nondimStateVars(&stateVarsCell[iQuadPt*_numVarsQuadPt],
			 _numVarsQuadPt);
      } // if

    } // for
    // Insert cell contribution into fields
    propertiesSection->updatePoint(*c_iter, &propertiesCell[0]);
    if (0 != _dbInitialState)
      stateVarsSection->updatePoint(*c_iter, &stateVarsCell[0]);
  } // for

  // Close databases
  _dbProperties->close();
  if (0 != _dbInitialState)
    _dbInitialState->close();
  logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
// Get properties
const pylith::topology::Field<pylith::topology::Mesh>&
pylith::materials::Material::getProperties() const
{ // getProperties
  return *_properties;
} // getProperties

// ----------------------------------------------------------------------
// Get state variables
const pylith::topology::Field<pylith::topology::Mesh>&
pylith::materials::Material::getStateVars() const
{ // getStateVars
  return *_stateVars;
} // getStateVars

// ----------------------------------------------------------------------
// Get physical property or state variable field.
void
pylith::materials::Material::getField(topology::Field<topology::Mesh> *field, const char* name) const
{ // getField
  assert(0 != field);
  assert(0 != _properties);
  assert(0 != _stateVars);

  int propertyIndex = -1;
  int stateVarIndex = -1;
  _findField(&propertyIndex, &stateVarIndex, name);

  // Get cell information
  const ALE::Obj<SieveMesh>& sieveMesh = field->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", _id);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  topology::FieldBase::VectorFieldEnum fieldType = topology::FieldBase::OTHER;
      
  if (propertyIndex >= 0) { // If field is a property
    int propOffset = 0;
    const string_vector& properties = _metadata.properties();
    assert(propertyIndex < properties.size());
    for (int i=0; i < propertyIndex; ++i)
      propOffset += 
	_metadata.fiberDim(properties[i].c_str(), Metadata::PROPERTY);
    const int fiberDim = _metadata.fiberDim(name, Metadata::PROPERTY);

    // Get properties section
    const ALE::Obj<RealSection>& propertiesSection = _properties->section();
    assert(!propertiesSection.isNull());
    const int totalPropsFiberDimLocal = (cells->size() > 0) ? 
      propertiesSection->getFiberDimension(*cells->begin()) : 0;
    int totalPropsFiberDim = 0;
    MPI_Allreduce((void *) &totalPropsFiberDimLocal, 
		  (void *) &totalPropsFiberDim, 1, 
		  MPI_INT, MPI_MAX, field->mesh().comm());
    assert(totalPropsFiberDim > 0);
    const int numPropsQuadPt = _numPropsQuadPt;
    const int numQuadPts = totalPropsFiberDim / numPropsQuadPt;
    assert(totalPropsFiberDim == numQuadPts * numPropsQuadPt);
    const int totalFiberDim = numQuadPts * fiberDim;

    // Get dimension scale information for properties.
    double_array propertyScales(numPropsQuadPt);
    propertyScales = 1.0;
    _dimProperties(&propertyScales[0], propertyScales.size());

    ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
    logger.stagePush("Materials");
    // Allocate buffer for property field if necessary.
    const ALE::Obj<RealSection>& fieldSection = field->section();
    bool useCurrentField = !fieldSection.isNull();
    if (!fieldSection.isNull()) {
      // check fiber dimension
      const int totalFiberDimCurrentLocal = (cells->size() > 0) ?
	fieldSection->getFiberDimension(*cells->begin()) : 0;
      int totalFiberDimCurrent = 0;
      MPI_Allreduce((void *) &totalFiberDimCurrentLocal, 
		    (void *) &totalFiberDimCurrent, 1, 
		    MPI_INT, MPI_MAX, field->mesh().comm());
      assert(totalFiberDimCurrent > 0);
      useCurrentField = totalFiberDim == totalFiberDimCurrent;
    } // if
    if (!useCurrentField) {
      field->newSection(cells, totalFiberDim);
      field->allocate();
    } // if
    assert(!fieldSection.isNull());
    field->label(name);
    field->scale(propertyScales[propOffset]);
    fieldType = _metadata.fieldType(name, Metadata::PROPERTY);
    logger.stagePop();
  
    // Buffer for property at cell's quadrature points
    double_array fieldCell(numQuadPts*fiberDim);
    double_array propertiesCell(numQuadPts*numPropsQuadPt);

    // Loop over cells
    for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
	 c_iter != cellsEnd;
	 ++c_iter) {
      propertiesSection->restrictPoint(*c_iter, 
				       &propertiesCell[0], propertiesCell.size());
   
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
	for (int i=0; i < fiberDim; ++i)
	  fieldCell[iQuad*fiberDim+i] = 
	    propertiesCell[iQuad*numPropsQuadPt+propOffset+i];

      fieldSection->updatePoint(*c_iter, &fieldCell[0]);
    } // for
  } else { // field is a state variable
    assert(stateVarIndex >= 0);
    
    int varOffset = 0;
    const string_vector& stateVars = _metadata.stateVars();
    assert(stateVarIndex < stateVars.size());
    for (int i=0; i < stateVarIndex; ++i)
      varOffset += 
	_metadata.fiberDim(stateVars[i].c_str(), Metadata::STATEVAR);
    const int fiberDim = _metadata.fiberDim(name, Metadata::STATEVAR);

    // Get state variables section
    const ALE::Obj<RealSection>& stateVarsSection = _stateVars->section();
    assert(!stateVarsSection.isNull());
    const int totalVarsFiberDimLocal = (cells->size() > 0) ?
      stateVarsSection->getFiberDimension(*cells->begin()) : 0;
    int totalVarsFiberDim = 0;
    MPI_Allreduce((void *) &totalVarsFiberDimLocal, 
		  (void *) &totalVarsFiberDim, 1, 
		  MPI_INT, MPI_MAX, field->mesh().comm());
    assert(totalVarsFiberDim > 0);
    const int numVarsQuadPt = _numVarsQuadPt;
    const int numQuadPts = totalVarsFiberDim / numVarsQuadPt;
    assert(totalVarsFiberDim == numQuadPts * numVarsQuadPt);
    const int totalFiberDim = numQuadPts * fiberDim;

    // Get dimension scale information for state variables.
    double_array stateVarScales(numVarsQuadPt);
    stateVarScales = 1.0;
    _dimStateVars(&stateVarScales[0], stateVarScales.size());

    ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
    logger.stagePush("Materials");
    // Allocate buffer for state variable field if necessary.
    const ALE::Obj<RealSection>& fieldSection = field->section();
    bool useCurrentField = !fieldSection.isNull();
    if (!fieldSection.isNull()) {
      // check fiber dimension
      const int totalFiberDimCurrentLocal = (cells->size() > 0) ?
	fieldSection->getFiberDimension(*cells->begin()) : 0;
      int totalFiberDimCurrent = 0;
      MPI_Allreduce((void *) &totalFiberDimCurrentLocal, 
		    (void *) &totalFiberDimCurrent, 1, 
		    MPI_INT, MPI_MAX, field->mesh().comm());
      assert(totalFiberDimCurrent > 0);
      useCurrentField = totalFiberDim == totalFiberDimCurrent;
    } // if
    if (!useCurrentField) {
      field->newSection(cells, totalFiberDim);
      field->allocate();
    } // if
    assert(!fieldSection.isNull());
    fieldType = _metadata.fieldType(name, Metadata::STATEVAR);
    field->label(name);
    field->scale(stateVarScales[varOffset]);
    logger.stagePop();
  
    // Buffer for state variable at cell's quadrature points
    double_array fieldCell(numQuadPts*fiberDim);
    double_array stateVarsCell(numQuadPts*numVarsQuadPt);
    
    // Loop over cells
    for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
	 c_iter != cellsEnd;
	 ++c_iter) {
      stateVarsSection->restrictPoint(*c_iter, 
				      &stateVarsCell[0], stateVarsCell.size());
      
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
	for (int i=0; i < fiberDim; ++i)
	  fieldCell[iQuad*fiberDim+i] = 
	    stateVarsCell[iQuad*numVarsQuadPt+varOffset+i];
      
      fieldSection->updatePoint(*c_iter, &fieldCell[0]);
    } // for
  } // if/else

  topology::FieldBase::VectorFieldEnum multiType = 
    topology::FieldBase::MULTI_OTHER;
  switch (fieldType)
    { // switch
    case topology::FieldBase::SCALAR:
      multiType = topology::FieldBase::MULTI_SCALAR;
      break;
    case topology::FieldBase::VECTOR:
      multiType = topology::FieldBase::MULTI_VECTOR;
      break;
    case topology::FieldBase::TENSOR:
      multiType = topology::FieldBase::MULTI_TENSOR;
      break;
    case topology::FieldBase::OTHER:
      multiType = topology::FieldBase::MULTI_OTHER;
      break;
    case topology::FieldBase::MULTI_SCALAR:
    case topology::FieldBase::MULTI_VECTOR:
    case topology::FieldBase::MULTI_TENSOR:
    case topology::FieldBase::MULTI_OTHER:
    default :
      std::cerr << "Bad vector field type for Material." << std::endl;
      assert(0);
      throw std::logic_error("Bad vector field type for Material.");
    } // switch
  field->vectorFieldType(multiType);
} // getField
  
// ----------------------------------------------------------------------
// Get indices for physical property or state variable field.
void
pylith::materials::Material::_findField(int* propertyIndex,
					int* stateVarIndex,
					const char* name) const
{ // _findField
  assert(0 != propertyIndex);
  assert(0 != stateVarIndex);

  *propertyIndex = -1;
  *stateVarIndex = -1;

  const std::string nameString = name;
  const string_vector& properties = _metadata.properties();
  const int numProperties = properties.size();
  for (int i=0; i < numProperties; ++i)
    if (nameString == properties[i]) {
      *propertyIndex = i;
      return;
    } // if

  const string_vector& stateVars = _metadata.stateVars();
  const int numStateVars = stateVars.size();
  for (int i=0; i < numStateVars; ++i)
    if (nameString == stateVars[i]) {
      *stateVarIndex = i;
      return;
    } // if

  if (propertyIndex < 0 && stateVarIndex < 0) {
    std::ostringstream msg;
    msg << "Unknown physical property or state variable '" << name
	<< "' for material '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // else
} // _findField
  

// End of file 
