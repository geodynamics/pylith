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

#include "BruneSlipFn.hh" // implementation of object methods

#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#include "pylith/utils/array.hh" // USES double_array

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::BruneSlipFn::BruneSlipFn(void) :
  _dbFinalSlip(0),
  _dbSlipTime(0),
  _dbPeakRate(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::BruneSlipFn::~BruneSlipFn(void)
{ // destructor
  _dbFinalSlip = 0;
  _dbSlipTime = 0;
  _dbPeakRate = 0;
} // destructor

// ----------------------------------------------------------------------
// Initialize slip time function.
void
pylith::faults::BruneSlipFn::initialize(const ALE::Obj<Mesh>& mesh,
					const ALE::Obj<Mesh>& faultMesh,
					const std::set<Mesh::point_type>& vertices,

					const spatialdata::geocoords::CoordSys* cs)
{ // initialize
  typedef std::set<Mesh::point_type>::const_iterator vert_iterator;  

  assert(!mesh.isNull());
  assert(!faultMesh.isNull());
  assert(0 != cs);
  assert(0 != _dbFinalSlip);
  assert(0 != _dbSlipTime);
  assert(0 != _dbPeakRate);

  const int spaceDim = cs->spaceDim();

  // Create and allocate sections for parameters
  delete _parameters; 
  _parameters = new topology::FieldsManager(mesh);
  if (0 == _parameters)
    throw std::runtime_error("Could not create manager for parameters of "
			     "Brune slip time function.");
  assert(0 != _parameters);
  
  // Parameter: final slip
  _parameters->addReal("final slip");
  const ALE::Obj<real_section_type>& finalSlip = 
    _parameters->getReal("final slip");
  assert(!finalSlip.isNull());

  // Parameter: slip initiation time
  _parameters->addReal("slip time");
  const ALE::Obj<real_section_type>& slipTime = 
    _parameters->getReal("slip time");
  assert(!slipTime.isNull());

  // Parameter: peak slip rate
  _parameters->addReal("peak rate");
  const ALE::Obj<real_section_type>& peakRate = 
    _parameters->getReal("peak rate");
  assert(!peakRate.isNull());

  const vert_iterator vBegin = vertices.begin();
  const vert_iterator vEnd = vertices.end();
  for (vert_iterator v_iter=vBegin; v_iter != vEnd; ++v_iter) {
    finalSlip->setFiberDimension(*v_iter, spaceDim);
    slipTime->setFiberDimension(*v_iter, 1);
    peakRate->setFiberDimension(*v_iter, 1);
  } // for
  mesh->allocate(finalSlip);
  mesh->allocate(slipTime);
  mesh->allocate(peakRate);
  
  // Open databases and set query values
  _dbFinalSlip->open();
  switch (spaceDim)
    { // switch
    case 1 : {
      const char* slipValues[] = {"slip"};
      _dbFinalSlip->queryVals(slipValues, 1);
      break;
    } // case 1
    case 2 : {
      const char* slipValues[] = {"slip", "fault-opening"};
      _dbFinalSlip->queryVals(slipValues, 2);
      break;
    } // case 2
    case 3 : {
      const char* slipValues[] = {"left-lateral-slip", "reverse-slip", 
				  "fault-opening"};
      _dbFinalSlip->queryVals(slipValues, 3);
      break;
    } // case 3
    default :
      assert(0);
    } // switch

  _dbSlipTime->open();
  const char* slipTimeValues[] = {"slip-time"};
  _dbSlipTime->queryVals(slipTimeValues, 1);

  _dbPeakRate->open();
  const char* peakRateValues[] = {"slip-rate"};
  _dbPeakRate->queryVals(peakRateValues, 1);

  // Get coordinates of vertices
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Query databases for parameters
  double_array slipData(spaceDim);
  double slipTimeData;
  double peakRateData;
  for (vert_iterator v_iter=vBegin; v_iter != vEnd; ++v_iter) {
    // Get coordinates of vertex
    const real_section_type::value_type* vCoords = 
      coordinates->restrictPoint(*v_iter);
    
    int err = _dbFinalSlip->query(&slipData[0], spaceDim, 
				  vCoords, spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find final slip at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoords[i];
      msg << ") using spatial database " << _dbFinalSlip->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    finalSlip->updatePoint(*v_iter, &slipData[0]);

    err = _dbSlipTime->query(&slipTimeData, 1, vCoords, spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find slip initiation time at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoords[i];
      msg << ") using spatial database " << _dbSlipTime->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    slipTime->updatePoint(*v_iter, &slipTimeData);

    err = _dbPeakRate->query(&peakRateData, 1, vCoords, spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find peak slip rate at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoords[i];
      msg << ") using spatial database " << _dbPeakRate->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    peakRate->updatePoint(*v_iter, &peakRateData);
  } // for

  // Close databases
  _dbFinalSlip->close();
  _dbSlipTime->close();
  _dbPeakRate->close();

  // Allocate slip field
  _slipField = new real_section_type(mesh->comm(), mesh->debug());
  for (vert_iterator v_iter=vBegin; v_iter != vEnd; ++v_iter)
    _slipField->setFiberDimension(*v_iter, spaceDim);
  mesh->allocate(_slipField);
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
const ALE::Obj<pylith::real_section_type>&
pylith::faults::BruneSlipFn::slip(const double t,
				  const std::set<Mesh::point_type>& vertices)
{ // slip
  typedef std::set<Mesh::point_type>::const_iterator vert_iterator;  

  assert(0 != _parameters);
  assert(!_slipField.isNull());
  
  // Get parameters
  const ALE::Obj<real_section_type>& finalSlip = 
    _parameters->getReal("final slip");
  assert(!finalSlip.isNull());

  const ALE::Obj<real_section_type>& slipTime = 
    _parameters->getReal("slip time");
  assert(!slipTime.isNull());

  const ALE::Obj<real_section_type>& peakRate = 
    _parameters->getReal("peak rate");
  assert(!peakRate.isNull());

  double_array slipValues(3);
  const vert_iterator vBegin = vertices.begin();
  const vert_iterator vEnd = vertices.end();
  for (vert_iterator v_iter=vBegin; v_iter != vEnd; ++v_iter) {
    // Get values of parameters at vertex
    const int numSlipValues = finalSlip->getFiberDimension(*v_iter);
    const real_section_type::value_type* vFinalSlip = 
      finalSlip->restrictPoint(*v_iter);
    const real_section_type::value_type* vSlipTime = 
      slipTime->restrictPoint(*v_iter);
    const real_section_type::value_type* vPeakRate = 
      peakRate->restrictPoint(*v_iter);

    double vFinalSlipMag = 0.0;
    for (int iSlip=0; iSlip < numSlipValues; ++iSlip)
      vFinalSlipMag += vFinalSlip[iSlip]*vFinalSlip[iSlip];
    vFinalSlipMag = sqrt(vFinalSlipMag);
    const double vSlip = _slip(t-vSlipTime[0], vFinalSlipMag, vPeakRate[0]);
    const double scale = vSlip / vFinalSlipMag;
    for (int iSlip=0; iSlip < numSlipValues; ++iSlip)
      slipValues[iSlip] = scale * vFinalSlip[iSlip];

    // Update field
    _slipField->updatePoint(*v_iter, &slipValues[0]);
  } // for

  return _slipField;
} // slip


// End of file 
