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

#include "pylith/feassemble/ParameterManager.hh" // USES ParameterManager
#include "pylith/feassemble/Quadrature.hh" // USES ParameterManager

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <petscmesh.h> // USES Mesh

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::Material::Material(void) :
  _parameters(0),
  _db(0),
  _id(0),
  _label("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::Material::~Material(void)
{ // destructor
  // Python db object owns database, so just set pointer to null
  _db = 0;

  delete _parameters; _parameters = 0;
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::materials::Material::Material(const Material& m) :
  _parameters(m._parameters),
  _db(m._db),
  _id(m._id),
  _label(m._label)
{ // copy constructor
} // copy constructor

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

  typedef ALE::Mesh::real_section_type real_section_type;
  typedef ALE::Mesh::topology_type topology_type;

  // Get cells associated with material
  const ALE::Mesh::int_section_type::patch_type patch = 0;
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<topology_type>& topology = coordinates->getTopology();
  const ALE::Obj<topology_type::label_sequence>& cells = 
    topology->getLabelStratum(patch, "material-id", _id);
  const topology_type::label_sequence::iterator cellsEnd = cells->end();

  // Check to make sure we have cells
  if (0 == cells->size()) {
    std::ostringstream msg;
    msg << "Could not find any cells for material '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if

  // Create sections to hold parameters for physical properties
  delete _parameters; _parameters = new feassemble::ParameterManager(mesh);
  const int numQuadPts = quadrature->numQuadPts();
  const int fiberDim = numQuadPts; // number of values in field per cell

  const int numParams = _numParameters();
  const char** paramNames = _parameterNames();

  ALE::Obj<real_section_type>* paramSections = 
    (numParams > 0) ? new ALE::Obj<real_section_type>[numParams] : 0;
  
  for (int iParam=0; iParam < numParams; ++iParam) {
    _parameters->addReal(paramNames[iParam]);
    paramSections[iParam] = _parameters->getReal(paramNames[iParam]);
    paramSections[iParam]->setFiberDimension(patch, cells, fiberDim);
    paramSections[iParam]->allocate();
  } // for

  // Setup database for querying
  const int numValues = _numDBValues();
  _db->open();
  _db->queryVals(_dbValues(), numValues);
  
  // Loop over cells
  double* queryData = (numValues > 0) ? new double[numValues] : 0;
  double* paramData = (numParams > 0) ? new double[numParams] : 0;
  double** cellData = (numParams > 0) ? new double*[numParams] : 0;
  for (int iParam = 0; iParam < numParams; ++iParam)
    cellData[iParam] = (numQuadPts > 0) ? new double[numQuadPts] : 0;
  for (topology_type::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    quadrature->computeGeometry(coordinates, *cellIter);

    const double* quadPts = quadrature->quadPts();
    const int spaceDim = quadrature->spaceDim();

    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, index=0; 
	 iQuadPt < numQuadPts; 
	 ++iQuadPt, index+=spaceDim) {
      const int err = _db->query(queryData, numValues, &quadPts[index],
				 spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find parameters for physical properties at \n"
	    << "(";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPts[index+spaceDim];
	msg << ") in material " << _label << "\n"
	    << "using spatial database " << _db->label() << ".";

	// Cleanup, then throw exception
	for (int iParam=0; iParam < numParams; ++iParam) {
	  delete[] cellData[iParam]; cellData[iParam] = 0;
	} // for
	delete[] cellData; cellData = 0;
	delete[] queryData; queryData = 0;
	delete[] paramData; paramData = 0;
	throw std::runtime_error(msg.str());
      } // if
      _dbToParameters(paramData, numParams, queryData, numValues);

      for (int iParam=0; iParam < numParams; ++iParam)
	cellData[iParam][iQuadPt] = paramData[iParam];
    } // for
    // Assemble cell contribution into fields
    for (int iParam=0; iParam < numParams; ++iParam)
      paramSections[iParam]->updateAdd(patch, *cellIter, cellData[iParam]);
  } // for
  for (int iParam=0; iParam < numParams; ++iParam) {
    delete[] cellData[iParam]; cellData[iParam] = 0;
  } // for
  delete[] cellData; cellData = 0;
  delete[] queryData; queryData = 0;
  delete[] paramData; paramData = 0;
  delete[] paramSections; paramSections = 0;

  // Close database
  _db->close();
} // initialize
  

// End of file 
