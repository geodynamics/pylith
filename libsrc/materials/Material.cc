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

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/utils/array.hh" // USES double_array, std::vector

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include "pylith/utils/sievetypes.hh" // USES Mesh

#include <string.h> // USES memcpy()
#include <strings.h> // USES strcasecmp()
#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::Material::Material(const int tensorSize,
				      const char** dbValues,
				      const char** initialStateDBValues,
				      const int numDBValues,
				      const PropMetaData* properties,
				      const int numProperties) :
  _dt(0.0),
  _totalPropsQuadPt(0),
  _dimension(0),
  _tensorSize(tensorSize),
  _initialStateSize(tensorSize),
  _initialStateDBValues(initialStateDBValues),
  _needNewJacobian(false),
  _db(0),
  _initialStateDB(0),
  _id(0),
  _label(""),
  _propMetaData(properties),
  _numProperties(numProperties),
  _dbValues(dbValues),
  _numDBValues(numDBValues)
{ // constructor
  for (int i=0; i < numProperties; ++i)
    _totalPropsQuadPt += properties[i].fiberDim;
  assert(_totalPropsQuadPt >= 0);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::Material::~Material(void)
{ // destructor
  // Python db object owns database, so just set pointer to null
  _db = 0;
  _initialStateDB = 0;
} // destructor

// ----------------------------------------------------------------------
// Get physical property parameters and initial state (if used) from database.
void
pylith::materials::Material::initialize(
			   const ALE::Obj<Mesh>& mesh,
			   const spatialdata::geocoords::CoordSys* cs,
			   pylith::feassemble::Quadrature* quadrature,
			   const spatialdata::units::Nondimensional& normalizer)
{ // initialize
  assert(0 != _db);
  assert(0 != cs);
  assert(0 != quadrature);
  assert(!mesh.isNull());

  // Get cells associated with material
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", _id);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsEnd = cells->end();

  // Create sections to hold physical properties and state variables.
  _properties = new real_section_type(mesh->comm(), mesh->debug());
  assert(!_properties.isNull());
  _properties->setChart(real_section_type::chart_type(*std::min_element(cells->begin(), cells->end()),
                                                      *std::max_element(cells->begin(), cells->end())+1));

  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();

  // Fiber dimension is number of quadrature points times number of
  // values per parameter
  const int totalPropsQuadPt = _totalPropsQuadPt;
  const int fiberDim = totalPropsQuadPt * numQuadPts;
  _properties->setFiberDimension(cells, fiberDim);
  mesh->allocate(_properties);

  // Fiber dimension for initial stresses is number of quadrature points times
  // initial state size.
  const int initialStateFiberDim = _initialStateSize * numQuadPts;

  // Container for data returned in query of initial state database
  const int initialStateSize = _initialStateSize;
  double_array initialStateQueryData(initialStateSize);

  // Container of initial state values at cell's quadrature points
  double_array initialStateCellData(initialStateFiberDim);

  // If initial state is being used, create a section to hold it.
  if (0 == _initialStateDB)
    assert(_initialState.isNull());
  else {
    _initialState = new real_section_type(mesh->comm(), mesh->debug());
    assert(!_initialState.isNull());
    _initialState->setChart(real_section_type::chart_type(*std::min_element(cells->begin(), cells->end()),
							   *std::max_element(cells->begin(), cells->end())+1));

    _initialState->setFiberDimension(cells, initialStateFiberDim);
    mesh->allocate(_initialState);

    // Setup database for querying
    _initialStateDB->open();
    _initialStateDB->queryVals(_initialStateDBValues, initialStateSize);
  } // if

  // Setup database for querying
  const int numValues = _numDBValues;
  _db->open();
  _db->queryVals(_dbValues, numValues);
  
  // Container for data returned in query of database
  double_array queryData(numValues);
  
  // Container of parameters at cell's quadrature points.
  double_array cellData(fiberDim);

  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    quadrature->computeGeometry(mesh, coordinates, *c_iter);

    const double_array& quadPts = quadrature->quadPts();

    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, index=0; 
	 iQuadPt < numQuadPts; 
	 ++iQuadPt, index+=spaceDim) {
      int err = _db->query(&queryData[0], numValues, &quadPts[index],
			   spaceDim, cs);

      if (err) {
	std::ostringstream msg;
	msg << "Could not find parameters for physical properties at \n"
	    << "(";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPts[index+i];
	msg << ") in material " << _label << "\n"
	    << "using spatial database '" << _db->label() << "'.";
	throw std::runtime_error(msg.str());
      } // if
      _dbToProperties(&cellData[totalPropsQuadPt*iQuadPt], queryData);
#if 0
      _nondimProperties(&cellData[totalPropsQuadPt*iQuadPt], totalPropsQuadPt);
#endif

      if (0 != _initialStateDB) {
	err = _initialStateDB->query(&initialStateQueryData[0], 
				     initialStateSize, 
				     &quadPts[index], spaceDim, cs);

	if (err) {
	  std::ostringstream msg;
	  msg << "Could not find initial state values at \n" << "(";
	  for (int i=0; i < spaceDim; ++i)
	    msg << "  " << quadPts[index+i];
	  msg << ") in material " << _label << "\n"
	      << "using spatial database '" << _initialStateDB->label() << "'.";
	  throw std::runtime_error(msg.str());
	} // if
#if 0 // nondimensionalize initial state
	_nondimInitialState(&initialStateCellData[initialStateSize*iQuadPt],
			    &initialStateQueryData[0], initialStateSize);
#else
	memcpy(&initialStateCellData[initialStateSize * iQuadPt],
	       &initialStateQueryData[0],
	       initialStateSize * sizeof(double));
#endif
      } // if

    } // for
    // Insert cell contribution into fields
    _properties->updatePoint(*c_iter, &cellData[0]);
    if (0 != _initialStateDB)
      _initialState->updatePoint(*c_iter, &initialStateCellData[0]);
  } // for

  // Close databases
  _db->close();
  if (0 != _initialStateDB)
    _initialStateDB->close();
} // initialize

// ----------------------------------------------------------------------
// Get type of field associated with physical property.
pylith::VectorFieldEnum
pylith::materials::Material::propertyFieldType(const char* name) const
{ // propertyFieldType
  VectorFieldEnum fieldType = OTHER_FIELD;

  // Find property in list of physical properties.
  int i=0;
  while (i < _numProperties)
    if (0 == strcasecmp(name, _propMetaData[i].name))
      break;
    else
      ++i;
  if (i < _numProperties)
    fieldType = _propMetaData[i].fieldType;
  else {
    std::ostringstream msg;
    msg << "Unknown physical property '" << name << "' for material '"
	<< _label << "'.";
    throw std::runtime_error(msg.str());
  } // else
 
  return fieldType;
} // propertyFieldType

// ----------------------------------------------------------------------
// Get physical property field.
void
pylith::materials::Material::propertyField(ALE::Obj<real_section_type>* field,
					   const char* name,
					   const ALE::Obj<Mesh>& mesh,
					   const int numQuadPts) const
{ // propertyField
  // Find property in list of physical properties.
  int i=0;
  int propOffset = 0;
  int fiberDim = 0;
  while (i < _numProperties)
    if (0 == strcasecmp(name, _propMetaData[i].name))
      break;
    else {
      propOffset += _propMetaData[i].fiberDim;
      ++i;
    } // else
  if (i < _numProperties)
    fiberDim = _propMetaData[i].fiberDim;
  else {
    std::ostringstream msg;
    msg << "Unknown physical property '" << name << "' for material '"
	<< _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  // Get cell information
  const ALE::Obj<Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", _id);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  
  const int totalFiberDim = numQuadPts * fiberDim;

  // Allocate buffer for property field.
  if (field->isNull() || 
      totalFiberDim != (*field)->getFiberDimension(*cells->begin())) {
    *field = new real_section_type(mesh->comm(), mesh->debug());
    (*field)->setChart(real_section_type::chart_type(*std::min_element(cells->begin(), cells->end()),
                                                     *std::max_element(cells->begin(), cells->end())+1));
    (*field)->setFiberDimension(cells, totalFiberDim);
    mesh->allocate(*field);
  } // if
  
  // Buffer for property at cell's quadrature points
  double_array fieldCell(fiberDim * numQuadPts);

  // Loop over cells
  const int totalPropsQuadPt = _totalPropsQuadPt;
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    const real_section_type::value_type* propVals = 
      _properties->restrictPoint(*c_iter);
   
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
#if 0
      _dimProperties(&propVals[totalPropsQuadPt*iQuad], totalPropsQuadPt);
#endif
      memcpy(&fieldCell[iQuad*fiberDim], 
	     &propVals[iQuad*totalPropsQuadPt+propOffset],
	     fiberDim*sizeof(double));
    } // for

    (*field)->updatePoint(*c_iter, &fieldCell[0]);
  } // for
} // propertyField
  

// End of file 
