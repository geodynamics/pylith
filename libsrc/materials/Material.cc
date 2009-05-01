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

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

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

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();

  // Get cells associated with material
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", _id);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);

  // Create field to hold physical properties.
  delete _properties; _properties = new topology::Field<topology::Mesh>(mesh);
  assert(0 != _properties);
  int fiberDim = numQuadPts * _numPropsQuadPt;
  _properties->newSection(cells, fiberDim);
  _properties->allocate();
  _properties->zero();
  const ALE::Obj<RealSection>& propertiesSection = _properties->section();
  assert(!propertiesSection.isNull());

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
    
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    quadrature->retrieveGeometry(*c_iter);

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
} // initialize

// ----------------------------------------------------------------------
// Get physical property or state variable field.
void
pylith::materials::Material::getField(topology::Field<topology::Mesh>* field,
				      const char* name) const
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
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  
  if (propertyIndex >= 0) { // If field is a property
    int propOffset = 0;
    const string_vector& properties = _metadata.properties();
    assert(propertyIndex < properties.size());
    for (int i=0; i < propertyIndex; ++i)
      propOffset += 
	_metadata.fiberDim(properties[i].c_str(), Metadata::PROPERTY);
    const int fiberDim = _metadata.fiberDim(name, Metadata::PROPERTY);

    // :TODO: Get scale information

    // Get properties section
    const ALE::Obj<RealSection>& propertiesSection = _properties->section();
    assert(!propertiesSection.isNull());
    const int totalPropsFiberDim = 
      propertiesSection->getFiberDimension(*cells->begin());
    const int numPropsQuadPt = _numPropsQuadPt;
    const int numQuadPts = totalPropsFiberDim / numPropsQuadPt;
    assert(totalPropsFiberDim == numQuadPts * numPropsQuadPt);
    const int totalFiberDim = numQuadPts * fiberDim;

    // Allocate buffer for property field.
    const ALE::Obj<RealSection>& fieldSection = field->section();
    if (fieldSection.isNull() ||
      totalFiberDim != fieldSection->getFiberDimension(*cells->begin())) {
      field->newSection(cells, totalFiberDim);
      field->allocate();
    } // if
    assert(!fieldSection.isNull());
    field->vectorFieldType(_metadata.fieldType(name, Metadata::PROPERTY));
    field->label(name);
  
    // Buffer for property at cell's quadrature points
    double_array fieldCell(numQuadPts*fiberDim);
    double_array propertiesCell(numQuadPts*numPropsQuadPt);

    // Loop over cells
    for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
	 c_iter != cellsEnd;
	 ++c_iter) {
      propertiesSection->restrictPoint(*c_iter, 
				       &propertiesCell[0], propertiesCell.size());
   
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	_dimProperties(&propertiesCell[iQuad*numPropsQuadPt], 
		       numPropsQuadPt);
	memcpy(&fieldCell[iQuad*fiberDim], 
	       &propertiesCell[iQuad*numPropsQuadPt+propOffset],
	       fiberDim*sizeof(double));
      } // for

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

    // :TODO: Get scale information

    // Get state variables section
    const ALE::Obj<RealSection>& stateVarsSection = _stateVars->section();
    assert(!stateVarsSection.isNull());
    const int totalVarsFiberDim = 
      stateVarsSection->getFiberDimension(*cells->begin());
    const int numVarsQuadPt = _numVarsQuadPt;
    const int numQuadPts = totalVarsFiberDim / numVarsQuadPt;
    assert(totalVarsFiberDim == numQuadPts * numVarsQuadPt);
    const int totalFiberDim = numQuadPts * fiberDim;

    // Allocate buffer for state variable field.
    const ALE::Obj<RealSection>& fieldSection = field->section();
    if (fieldSection.isNull() ||
	totalFiberDim != fieldSection->getFiberDimension(*cells->begin())) {
      field->newSection(cells, totalFiberDim);
      field->allocate();
    } // if
    assert(!fieldSection.isNull());
    field->vectorFieldType(_metadata.fieldType(name, Metadata::STATEVAR));
    field->label(name);
  
    // Buffer for state variable at cell's quadrature points
    double_array fieldCell(numQuadPts*fiberDim);
    double_array stateVarsCell(numQuadPts*numVarsQuadPt);
    
    // Loop over cells
    for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
	 c_iter != cellsEnd;
	 ++c_iter) {
      stateVarsSection->restrictPoint(*c_iter, 
				      &stateVarsCell[0], stateVarsCell.size());
      
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	_dimStateVars(&stateVarsCell[iQuad*numVarsQuadPt], 
		      numVarsQuadPt);
	memcpy(&fieldCell[iQuad*fiberDim], 
	       &stateVarsCell[iQuad*numVarsQuadPt+varOffset],
	       fiberDim*sizeof(double));
      } // for

      fieldSection->updatePoint(*c_iter, &fieldCell[0]);
    } // for
  } // if/else
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
