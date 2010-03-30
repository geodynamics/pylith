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

#include "FrictionModel.hh" // implementation of object methods

#include "pylith/topology/SubMesh.hh" // USES Mesh
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
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::RestrictVisitor RestrictVisitor;

// ----------------------------------------------------------------------
// Default constructor.
pylith::friction::FrictionModel::FrictionModel(const materials::Metadata& metadata) :
  _dt(0.0),
  _properties(0),
  _stateVars(0),
  _normalizer(new spatialdata::units::Nondimensional),
  _numPropsVertex(0),
  _numVarsVertex(0),
  _dbProperties(0),
  _dbInitialState(0),
  _label(""),
  _metadata(metadata)
{ // constructor
  const string_vector& properties = metadata.properties();
  const int numProperties = properties.size();
  for (int i=0; i < numProperties; ++i)
    _numPropsVertex += metadata.fiberDim(properties[i].c_str(),
					 materials::Metadata::PROPERTY);
  assert(_numPropsVertex >= 0);

  const string_vector& stateVars = metadata.stateVars();
  const int numStateVars = stateVars.size();
  for (int i=0; i < numStateVars; ++i)
    _numVarsVertex += metadata.fiberDim(stateVars[i].c_str(),
				  materials::Metadata::STATEVAR);
  assert(_numVarsVertex >= 0);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::friction::FrictionModel::~FrictionModel(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::friction::FrictionModel::deallocate(void)
{ // deallocate
  delete _normalizer; _normalizer = 0;
  delete _properties; _properties = 0;
  delete _stateVars; _stateVars = 0;

  _dbProperties = 0; // :TODO: Use shared pointer.
  _dbInitialState = 0; // :TODO: Use shared pointer.
} // deallocate
  
// ----------------------------------------------------------------------
// Set scales used to nondimensionalize physical properties.
void
pylith::friction::FrictionModel::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
  if (0 == _normalizer)
    _normalizer = new spatialdata::units::Nondimensional(dim);
  else
    *_normalizer = dim;
} // normalizer

// ----------------------------------------------------------------------
// Get physical property parameters and initial state (if used) from database.
void
pylith::friction::FrictionModel::initialize(
			const topology::SubMesh& faultMesh,
			feassemble::Quadrature<topology::SubMesh>* quadrature,
			const topology::Field<topology::SubMesh>& area)
{ // initialize
  assert(0 != _dbProperties);
  assert(0 != quadrature);

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //logger.stagePush("Friction");

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int numBasis = quadrature->numBasis();
  const int spaceDim = quadrature->spaceDim();
  const double_array& quadWts = quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  double_array quadPtsGlobal(numQuadPts*spaceDim);

  // Get cells associated with friction interface
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh.sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
    faultSieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin =
    vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(0 != cs);

#if !defined(PRECOMPUTE_GEOMETRY)
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    faultSieveMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
        coordinatesCell.size(), &coordinatesCell[0]);
#endif

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();

  // Query database for properties

  // Create arrays for querying.
  const int numDBProperties = _metadata.numDBProperties();
  double_array propertiesDBQuery(numDBProperties);
  double_array propertiesQuadPt(_numPropsVertex);

  // Create field to hold physical properties.
  delete _properties; _properties = 
			new topology::Field<topology::SubMesh>(faultMesh);
  _properties->label("properties");
  assert(0 != _properties);
  int fiberDim = _numPropsVertex;
  _properties->newSection(vertices, fiberDim);
  _properties->allocate();
  _properties->zero();
  const ALE::Obj<RealSection>& propertiesSection = _properties->section();
  assert(!propertiesSection.isNull());
  double_array propertiesCell(numBasis*numDBProperties);
  topology::Mesh::UpdateAddVisitor propertiesVisitor(*propertiesSection,
        &propertiesCell[0]);

  // Setup database for querying for physical properties
  assert(0 != _dbProperties);
  _dbProperties->open();
  _dbProperties->queryVals(_metadata.dbProperties(),
			   _metadata.numDBProperties());

  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    propertiesCell = 0.0;

    const double_array& quadPtsNonDim = quadrature->quadPts();
    quadPtsGlobal = quadPtsNonDim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
        lengthScale);

    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, index=0;
        iQuadPt < numQuadPts;
        ++iQuadPt, index+=spaceDim) {
      int err = _dbProperties->query(&propertiesDBQuery[0], numDBProperties,
          &quadPtsGlobal[index], spaceDim, cs);
      if (err) {
        std::ostringstream msg;
        msg << "Could not find parameters for physical properties at \n" << "(";
        for (int i = 0; i < spaceDim; ++i)
          msg << "  " << quadPtsGlobal[index + i];
        msg << ") in friction model " << _label << "\n"
            << "using spatial database '" << _dbProperties->label() << "'.";
        throw std::runtime_error(msg.str());
      } // if
      _dbToProperties(&propertiesQuadPt[0], propertiesDBQuery);

      // Get cell geometry information that depends on cell
      const double_array& basis = quadrature->basis();
      const double_array& jacobianDet = quadrature->jacobianDet();

      // Compute properties weighted by area
      const double wt = quadWts[iQuadPt] * jacobianDet[iQuadPt];
      for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
        const double dArea = wt * basis[iQuadPt*numBasis+iBasis];
        for (int iProp = 0; iProp < numDBProperties; ++iProp)
          propertiesCell[iBasis*numDBProperties+iProp]
              += propertiesQuadPt[iProp] * dArea;
      } // for
    } // for
    propertiesVisitor.clear();
    faultSieveMesh->updateClosure(*c_iter, propertiesVisitor);
  } // for
  // Close properties database
  _dbProperties->close();

  _properties->complete(); // Assemble contributions

  // Loop over vertices and divide by area to get weighted values and
  // nondimensionalize properties.
  const ALE::Obj<RealSection>& areaSection = area.section();
  assert(!areaSection.isNull());

  double_array propertiesVertex(_numPropsVertex);
  double areaVertex = 0.0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    propertiesSection->restrictPoint(*v_iter,
        &propertiesVertex[0], propertiesVertex.size());
    _nondimProperties(&propertiesVertex[0], _numPropsVertex);
    areaSection->restrictPoint(*v_iter, &areaVertex, 1);
    assert(areaVertex > 0.0);
    propertiesVertex /= areaVertex;
    propertiesSection->updatePoint(*v_iter, &propertiesVertex[0]);
  } // for

  // Create field to hold state variables. We create the field even
  // if there is no initial state, because this we will use this field
  // to hold the state variables.
  delete _stateVars; _stateVars = new topology::Field<topology::SubMesh>(faultMesh);
  _stateVars->label("state variables");
  fiberDim = _numVarsVertex;
  if (fiberDim > 0) {
    assert(0 != _stateVars);
    assert(0 != _properties);
    _stateVars->newSection(*_properties, fiberDim);
    _stateVars->allocate();
    _stateVars->zero();
  } // if

  // Query database for initial state variables
  if (0 != _dbInitialState) {
    assert(_numVarsVertex > 0);

    // Create arrays for querying
    const int numDBStateVars = _metadata.numDBStateVars();
    double_array stateVarsDBQuery(numDBStateVars);
    double_array stateVarsQuadPt(_numVarsVertex);

    // Setup database for querying for initial state variables
    _dbInitialState->open();
    _dbInitialState->queryVals(_metadata.dbStateVars(),
             _metadata.numDBStateVars());

    const ALE::Obj<RealSection>& stateVarsSection =_stateVars->section();
    assert(!stateVarsSection.isNull());
    double_array stateVarsCell(numBasis*numDBStateVars);
    topology::Mesh::UpdateAddVisitor stateVarsVisitor(*stateVarsSection,
                                                      &stateVarsCell[0]);

    for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
        c_iter != cellsEnd;
        ++c_iter) {
      // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
      quadrature->retrieveGeometry(*c_iter);
#else
      coordsVisitor.clear();
      faultSieveMesh->restrictClosure(*c_iter, coordsVisitor);
      quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

      const double_array& quadPtsNonDim = quadrature->quadPts();
      quadPtsGlobal = quadPtsNonDim;
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
          lengthScale);

      stateVarsCell = 0.0;

      // Loop over quadrature points in cell and query database
      for (int iQuadPt=0, index=0;
          iQuadPt < numQuadPts;
          ++iQuadPt, index+=spaceDim) {
        int err = _dbInitialState->query(&stateVarsDBQuery[0], numDBStateVars,
          &quadPtsGlobal[index], spaceDim, cs);
        if (err) {
          std::ostringstream msg;
          msg << "Could not find initial state variables at \n" << "(";
          for (int i = 0; i < spaceDim; ++i)
            msg << "  " << quadPtsGlobal[index + i];
          msg << ") in friction model " << _label << "\n"
              << "using spatial database '" << _dbInitialState->label() << "'.";
          throw std::runtime_error(msg.str());
        } // if
        _dbToStateVars(&stateVarsQuadPt[0], stateVarsDBQuery);

        // Get cell geometry information that depends on cell
        const double_array& basis = quadrature->basis();
        const double_array& jacobianDet = quadrature->jacobianDet();

        // Compute state variables weighted by area
        const double wt = quadWts[iQuadPt] * jacobianDet[iQuadPt];
        for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
          const double dArea = wt * basis[iQuadPt*numBasis+iBasis];
          for (int iVar = 0; iVar < numDBStateVars; ++iVar)
            stateVarsCell[iBasis*numDBStateVars+iVar]
                += stateVarsDBQuery[iVar] * dArea;
        } // for
      } // for
      stateVarsVisitor.clear();
      faultSieveMesh->updateClosure(*c_iter, stateVarsVisitor);
    } // for
    // Close database
    _dbInitialState->close();

    _stateVars->complete(); // Assemble contributions.

    // Loop over vertices and divide by area to get weighted values and
    // nondimensionalize properties.
    double_array stateVarsVertex(_numVarsVertex);

    for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
         v_iter != verticesEnd;
         ++v_iter) {
      stateVarsSection->restrictPoint(*v_iter,
          &stateVarsVertex[0], stateVarsVertex.size());
      _nondimStateVars(&stateVarsVertex[0], _numVarsVertex);
      areaSection->restrictPoint(*v_iter, &areaVertex, 1);
      assert(areaVertex > 0.0);
      stateVarsVertex /= areaVertex;
      stateVarsSection->updatePoint(*v_iter, &stateVarsVertex[0]);
    } // for
  } // if

  // Setup buffers for restrict/update of properties and state variables.
  _propertiesVertex.resize(_numPropsVertex);
  _stateVarsVertex.resize(_numVarsVertex);

  //logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
// Get the properties field.
const pylith::topology::Field<pylith::topology::SubMesh>*
pylith::friction::FrictionModel::propertiesField() const
{ // propertiesField
  return _properties;
} // propertiesField

// ----------------------------------------------------------------------
// Get the state variables field.
const pylith::topology::Field<pylith::topology::SubMesh>*
pylith::friction::FrictionModel::stateVarsField() const
{ // stateVarsField
  return _stateVars;
} // stateVarsField

// ----------------------------------------------------------------------
// Check whether material has a field as a property.
bool
pylith::friction::FrictionModel::hasProperty(const char* name)
{ // hasProperty
  int propertyIndex = -1;
  int stateVarIndex = -1;
  _findField(&propertyIndex, &stateVarIndex, name);
  return (propertyIndex >= 0);
} // hasProperty

// ----------------------------------------------------------------------
// Check whether material has a field as a state variable.
bool
pylith::friction::FrictionModel::hasStateVar(const char* name)
{ // hasStateVar
  int propertyIndex = -1;
  int stateVarIndex = -1;
  _findField(&propertyIndex, &stateVarIndex, name);
  return (stateVarIndex >= 0);
} // hasStateVar

// ----------------------------------------------------------------------
// Get physical property or state variable field.
void
pylith::friction::FrictionModel::getField(topology::Field<topology::SubMesh> *field,
					  const char* name) const
{ // getField
  // Logging of allocation is handled by getField() caller since it
  // manages the memory for the field argument.

  assert(0 != field);
  assert(0 != _properties);
  assert(0 != _stateVars);

  int propertyIndex = -1;
  int stateVarIndex = -1;
  _findField(&propertyIndex, &stateVarIndex, name);
  if (propertyIndex < 0 && stateVarIndex < 0) {
    std::ostringstream msg;
    msg << "Unknown physical property or state variable '" << name
	<< "' for material '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  // Get cell information
  const ALE::Obj<SieveSubMesh>& sieveSubMesh = field->mesh().sieveMesh();
  assert(!sieveSubMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
    sieveSubMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin =
    vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

  topology::FieldBase::VectorFieldEnum fieldType = topology::FieldBase::OTHER;
      
  if (propertyIndex >= 0) { // If field is a property
    int propOffset = 0;
    const string_vector& properties = _metadata.properties();
    assert(propertyIndex < properties.size());
    for (int i = 0; i < propertyIndex; ++i)
      propOffset += _metadata.fiberDim(properties[i].c_str(),
          materials::Metadata::PROPERTY);
    const int fiberDim =
        _metadata.fiberDim(name, materials::Metadata::PROPERTY);

    // Get properties section
    const ALE::Obj<RealSection>& propertiesSection = _properties->section();
    assert(!propertiesSection.isNull());
    const int numPropsVertex = _numPropsVertex;

    // Get dimension scale information for properties.
    double_array propertyScales(numPropsVertex);
    propertyScales = 1.0;
    _dimProperties(&propertyScales[0], propertyScales.size());

    // Allocate buffer for property field if necessary.
    const ALE::Obj<RealSection>& fieldSection = field->section();
    bool useCurrentField = !fieldSection.isNull();
    if (!fieldSection.isNull()) {
      // check fiber dimension
      const int
          totalFiberDimCurrentLocal =
              (vertices->size() > 0) ? fieldSection->getFiberDimension(
                  *verticesBegin) : 0;
      int totalFiberDimCurrent = 0;
      MPI_Allreduce((void *) &totalFiberDimCurrentLocal,
          (void *) &totalFiberDimCurrent, 1,
          MPI_INT, MPI_MAX, field->mesh().comm());
      assert(totalFiberDimCurrent > 0);
      useCurrentField = fiberDim == totalFiberDimCurrent;
    } // if
    if (!useCurrentField) {
      ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
      logger.stagePush("Output");
      field->newSection(vertices, fiberDim);
      field->allocate();
      logger.stagePop();
    } // if
    assert(!fieldSection.isNull());
    field->label(name);
    field->scale(propertyScales[propOffset]);
    fieldType = _metadata.fieldType(name, materials::Metadata::PROPERTY);

    // Buffer for property at cell's quadrature points
    double_array fieldVertex(fiberDim);
    double_array propertiesVertex(numPropsVertex);

    // Loop over cells
    for (SieveSubMesh::label_sequence::iterator v_iter = verticesBegin;
        v_iter != verticesEnd;
        ++v_iter) {
      propertiesSection->restrictPoint(*v_iter, &propertiesVertex[0],
          propertiesVertex.size());

      for (int i = 0; i < fiberDim; ++i)
        fieldVertex[i] = propertiesVertex[propOffset+i];

      fieldSection->updatePoint(*v_iter, &fieldVertex[0]);
    } // for
  } else { // field is a state variable
    assert(stateVarIndex >= 0);

    int varOffset = 0;
    const string_vector& stateVars = _metadata.stateVars();
    assert(stateVarIndex < stateVars.size());
    for (int i = 0; i < stateVarIndex; ++i)
      varOffset += _metadata.fiberDim(stateVars[i].c_str(),
          materials::Metadata::STATEVAR);
    const int fiberDim =
        _metadata.fiberDim(name, materials::Metadata::STATEVAR);

    // Get state variables section
    const ALE::Obj<RealSection>& stateVarsSection = _stateVars->section();
    assert(!stateVarsSection.isNull());
    const int numVarsVertext = _numVarsVertex;

    // Get dimension scale information for state variables.
    double_array stateVarScales(_numVarsVertex);
    stateVarScales = 1.0;
    _dimStateVars(&stateVarScales[0], stateVarScales.size());

    // Allocate buffer for state variable field if necessary.
    const ALE::Obj<RealSection>& fieldSection = field->section();
    bool useCurrentField = !fieldSection.isNull();
    if (!fieldSection.isNull()) {
      // check fiber dimension
      const int
          totalFiberDimCurrentLocal =
              (vertices->size() > 0) ? fieldSection->getFiberDimension(
                  *verticesBegin) : 0;
      int totalFiberDimCurrent = 0;
      MPI_Allreduce((void *) &totalFiberDimCurrentLocal,
          (void *) &totalFiberDimCurrent, 1,
          MPI_INT, MPI_MAX, field->mesh().comm());
      assert(totalFiberDimCurrent > 0);
      useCurrentField = fiberDim == totalFiberDimCurrent;
    } // if
    if (!useCurrentField) {
      ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
      logger.stagePush("Output");
      field->newSection(vertices, fiberDim);
      field->allocate();
      logger.stagePop();
    } // if
    assert(!fieldSection.isNull());
    fieldType = _metadata.fieldType(name, materials::Metadata::STATEVAR);
    field->label(name);
    field->scale(stateVarScales[varOffset]);

    // Buffer for state variable at cell's quadrature points
    double_array fieldVertex(fiberDim);
    double_array stateVarsVertex(_numVarsVertex);

    // Loop over cells
    for (SieveSubMesh::label_sequence::iterator v_iter = verticesBegin;
          v_iter != verticesEnd;
          ++v_iter) {
      stateVarsSection->restrictPoint(*v_iter, &stateVarsVertex[0],
          stateVarsVertex.size());

      for (int i = 0; i < fiberDim; ++i)
        fieldVertex[i] = stateVarsVertex[varOffset + i];

      fieldSection->updatePoint(*v_iter, &fieldVertex[0]);
    } // for
  } // if/else

  topology::FieldBase::VectorFieldEnum multiType =
      topology::FieldBase::MULTI_OTHER;
  switch (fieldType) { // switch
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
  default:
    std::cerr << "Bad vector field type '" << fieldType << "'." << std::endl;
    assert(0);
    throw std::logic_error("Bad vector field type for FrictionModel.");
  } // switch
  field->vectorFieldType(multiType);
} // getField
  
// ----------------------------------------------------------------------
// Retrieve parameters for physical properties and state variables
// for vertex.
void
pylith::friction::FrictionModel::retrievePropsAndVars(const int vertex)
{ // retrievePropsAndVars
  assert(0 != _properties);
  assert(0 != _stateVars);

  const ALE::Obj<RealSection>& propertiesSection = _properties->section();
  assert(!propertiesSection.isNull());
  assert(_propertiesVertex.size() == propertiesSection->getFiberDimension(vertex));
  propertiesSection->restrictPoint(vertex, &_propertiesVertex[0],
           _propertiesVertex.size());

  if (_numVarsVertex > 0) {
    const ALE::Obj<RealSection>& stateVarsSection = _stateVars->section();
    assert(!stateVarsSection.isNull());
    assert(_stateVarsVertex.size() == stateVarsSection->getFiberDimension(vertex));
    stateVarsSection->restrictPoint(vertex, &_stateVarsVertex[0],
            _stateVarsVertex.size());
  } // if
} // retrievePropsAndVars

// ----------------------------------------------------------------------
// Compute friction at vertex.
double
pylith::friction::FrictionModel::calcFriction(const double slip,
                                              const double slipRate,
                                              const double normalTraction)
{ // calcFriction
  const double friction =
    _calcFriction(slip, slipRate, normalTraction,
		  &_propertiesVertex[0], _propertiesVertex.size(),
		  &_stateVarsVertex[0], _stateVarsVertex.size());

  return friction;
} // calcFriction

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::FrictionModel::updateStateVars(const double slip,
						 const double slipRate,
						 const double normalTraction)
{ // updateStateVars
  _updateStateVars(slip, slipRate, normalTraction,
		   &_stateVarsVertex[0], _stateVarsVertex.size(),
		   &_propertiesVertex[0], _propertiesVertex.size());
} // updateStateVars

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::FrictionModel::_updateStateVars(const double slip,
    const double slipRate,
    const double normalTraction,
    double* const stateVars,
    const int numStateVars,
    const double* properties,
    const int numProperties)
{ // _updateStateVars
} // _updateStateVars

// ----------------------------------------------------------------------
// Get indices for physical property or state variable field.
void
pylith::friction::FrictionModel::_findField(int* propertyIndex,
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
} // _findField
  

// End of file 
