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

#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/utils/array.hh" // USES double_array, std::vector

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include "pylith/utils/sievetypes.hh" // USES Mesh

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::Material::Material(const int* numParamValues,
				      const int size) :
  _dt(0.0),
  _numParamsQuadPt(0),
  _dimension(0),
  _needNewJacobian(false),
  _db(0),
  _id(0),
  _label(""),
  _numParamValues(numParamValues, size)
{ // constructor
  for (int i=0; i < size; ++i)
    _numParamsQuadPt += numParamValues[i];
  assert(_numParamsQuadPt >= 0);
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
  _parameters = new real_section_type(mesh->comm(), mesh->debug());
  assert(!_parameters.isNull());
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();

  const int numParams = _numParamValues.size();

  // Fiber dimension is number of quadrature points times number of
  // values per parameter
  const int numParamsQuadPt = _numParamsQuadPt;
  const int fiberDim = numParamsQuadPt * numQuadPts;
  _parameters->setFiberDimension(cells, fiberDim);
  mesh->allocate(_parameters);

  // Setup database for querying
  const int numValues = _numDBValues();
  _db->open();
  _db->queryVals(_dbValues(), numValues);
  
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
      _dbToParameters(&cellData[numParamsQuadPt*iQuadPt], numParamsQuadPt, 
		      queryData);

    } // for
    // Insert cell contribution into fields
    _parameters->updatePoint(*c_iter, &cellData[0]);
  } // for

  // Close database
  _db->close();
} // initialize
  

// End of file 
