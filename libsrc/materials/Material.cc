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
#include "pylith/utils/array.hh" // USES double_array, std::vector

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <petscmesh.h> // USES Mesh

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::Material::Material(void) :
  _parameters(0),
  _dimension(0),
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
  _dimension(m._dimension),
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
  assert(!mesh.isNull());

  // Get cells associated with material
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<ALE::Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", _id);
  assert(!cells.isNull());
  const ALE::Mesh::label_sequence::iterator cellsEnd = cells->end();

  // Check to make sure we have cells
  if (0 == cells->size()) {
    std::ostringstream msg;
    msg << "Could not find any cells for material '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if

  // Create sections to hold parameters for physical properties
  delete _parameters; _parameters = new feassemble::ParameterManager(mesh);
  assert(0 != _parameters);
  const int numQuadPts = quadrature->numQuadPts();
  const int fiberDim = numQuadPts; // number of values in field per cell

  const int numParams = _numParameters();
  const char** paramNames = _parameterNames();

  std::vector<ALE::Obj<real_section_type> > paramSections(numParams);
  
  for (int iParam=0; iParam < numParams; ++iParam) {
    _parameters->addReal(paramNames[iParam]);
    paramSections[iParam] = _parameters->getReal(paramNames[iParam]);
    assert(!paramSections[iParam].isNull());
    paramSections[iParam]->setFiberDimension(cells, fiberDim);
    mesh->allocate(paramSections[iParam]);
  } // for

  // Setup database for querying
  const int numValues = _numDBValues();
  _db->open();
  _db->queryVals(_dbValues(), numValues);
  
  // Loop over cells
  double_array queryData(numValues);
  double_array paramData(numParams);
  std::vector<double_array> cellData(numParams);
  for (int iParam = 0; iParam < numParams; ++iParam)
    cellData[iParam].resize(numQuadPts);
  for (ALE::Mesh::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    quadrature->computeGeometry(mesh, coordinates, *cellIter);

    const double_array& quadPts = quadrature->quadPts();
    const int spaceDim = quadrature->spaceDim();

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
	  msg << "  " << quadPts[index+spaceDim];
	msg << ") in material " << _label << "\n"
	    << "using spatial database " << _db->label() << ".";
	throw std::runtime_error(msg.str());
      } // if
      _dbToParameters(&paramData, queryData);

      for (int iParam=0; iParam < numParams; ++iParam)
	cellData[iParam][iQuadPt] = paramData[iParam];
    } // for
    // Assemble cell contribution into fields
    for (int iParam=0; iParam < numParams; ++iParam)
      mesh->updateAdd(paramSections[iParam], *cellIter, &cellData[iParam][0]);
  } // for

  // Close database
  _db->close();
} // initialize
  

// End of file 
