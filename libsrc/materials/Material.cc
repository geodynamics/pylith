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

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::Material::Material(const char** dbValues,
				      const int numDBValues,
				      const PropMetaData* properties,
				      const int numProperties) :
  _dt(0.0),
  _totalPropsQuadPt(0),
  _dimension(0),
  _needNewJacobian(false),
  _db(0),
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
} // destructor

// ----------------------------------------------------------------------
// Get physical property parameters from database.
void
pylith::materials::Material::initialize(const ALE::Obj<ALE::Mesh>& mesh,
					const spatialdata::geocoords::CoordSys* cs,
					pylith::feassemble::Quadrature* quadrature)
{ // initialize
  assert(0 != _db);
  assert(0 != cs);
  assert(0 != quadrature);
  assert(!mesh.isNull());

  // Get cells associated with material
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<ALE::Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", _id);
  assert(!cells.isNull());
  const ALE::Mesh::label_sequence::iterator cellsEnd = cells->end();

  // Create sections to hold parameters for physical properties
  _properties = new real_section_type(mesh->comm(), mesh->debug());
  assert(!_properties.isNull());
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();

  // Fiber dimension is number of quadrature points times number of
  // values per parameter
  const int totalPropsQuadPt = _totalPropsQuadPt;
  const int fiberDim = totalPropsQuadPt * numQuadPts;
  _properties->setFiberDimension(cells, fiberDim);
  mesh->allocate(_properties);

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
      const int err = _db->query(&queryData[0], numValues, &quadPts[index],
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

    } // for
    // Insert cell contribution into fields
    _properties->updatePoint(*c_iter, &cellData[0]);
  } // for

  // Close database
  _db->close();
} // initialize

// ----------------------------------------------------------------------
// Get metadata for physical property. Values are returned through the
// arguments.
void
pylith::materials::Material::propertyInfo(int* space,
					  int* fiberDim,
					  VectorFieldEnum* fieldType,
					  const char* name) const
{ // propertyInfo
  int i=0;
  while (i < _numProperties)
    if (0 == strcasecmp(name, _propMetaData[i].name))
      break;
  if (i < _numProperties) {
    *space = i;
    *fiberDim = _propMetaData[i].fiberDim;
    *fieldType = _propMetaData[i].fieldType;
  } else {
    std::ostringstream msg;
    msg << "Unknown physical property '" << name << "' for material '"
	<< _label << "'.";
    throw std::runtime_error(msg.str());
  } // else
} // propertyInfo
  

// End of file 
