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

namespace pylith {
  namespace faults {
    namespace _BruneSlipFn {
      const int offsetPeakRate = 0;
      const int offsetSlipTime = 1;
    } // _BruneSlipFn
  } // faults
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::BruneSlipFn::BruneSlipFn(void) :
  _dbFinalSlip(0),
  _dbSlipTime(0),
  _dbPeakRate(0),
  _spaceDim(0)
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
pylith::faults::BruneSlipFn::initialize(
				 const ALE::Obj<Mesh>& faultMesh,
				 const spatialdata::geocoords::CoordSys* cs)
{ // initialize
  assert(!faultMesh.isNull());
  assert(0 != cs);
  assert(0 != _dbFinalSlip);
  assert(0 != _dbSlipTime);
  assert(0 != _dbPeakRate);

  _spaceDim = cs->spaceDim();
  const int spaceDim = _spaceDim;
  const int indexFinalSlip = 0;
  const int indexPeakRate = spaceDim + _BruneSlipFn::offsetPeakRate;
  const int indexSlipTime = spaceDim + _BruneSlipFn::offsetSlipTime;

  // Get vertices in fault mesh
  const ALE::Obj<Mesh::label_sequence>& vertices = faultMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  const int fiberDim = spaceDim + 2;
  _parameters = new real_section_type(faultMesh->comm(), faultMesh->debug());
  _parameters->addSpace(); // final slip
  _parameters->addSpace(); // peak slip rate
  _parameters->addSpace(); // slip time
  assert(3 == _parameters->getNumSpaces());
  _parameters->setFiberDimension(vertices, fiberDim);
  _parameters->setFiberDimension(vertices, spaceDim, 0); // final slip
  _parameters->setFiberDimension(vertices, 1, 1); // peak slip rate
  _parameters->setFiberDimension(vertices, 1, 2); // slip time
  faultMesh->allocate(_parameters);
  assert(!_parameters.isNull());

  // Open databases and set query values
  _dbFinalSlip->open();
  switch (spaceDim)
    { // switch
    case 1 : {
      const char* slipValues[] = {"fault-opening"};
      _dbFinalSlip->queryVals(slipValues, 1);
      break;
    } // case 1
    case 2 : {
      const char* slipValues[] = {"left-lateral-slip", "fault-opening"};
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
    faultMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  double_array paramsVertex(fiberDim);

  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter) {

    // Get coordinates of vertex
    const real_section_type::value_type* coordsVertex = 
      coordinates->restrictPoint(*v_iter);
    assert(0 != coordsVertex);
    
    int err = _dbFinalSlip->query(&paramsVertex[indexFinalSlip], spaceDim, 
				  coordsVertex, spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find final slip at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << coordsVertex[i];
      msg << ") using spatial database " << _dbFinalSlip->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    err = _dbPeakRate->query(&paramsVertex[indexPeakRate], 1, 
			     coordsVertex, spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find peak slip rate at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << coordsVertex[i];
      msg << ") using spatial database " << _dbPeakRate->label() << ".";
      throw std::runtime_error(msg.str());
    } // if

    err = _dbSlipTime->query(&paramsVertex[indexSlipTime], 1, 
			     coordsVertex, spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find slip initiation time at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << coordsVertex[i];
      msg << ") using spatial database " << _dbSlipTime->label() << ".";
      throw std::runtime_error(msg.str());
    } // if

    _parameters->updatePoint(*v_iter, &paramsVertex[0]);
  } // for

  // Close databases
  _dbFinalSlip->close();
  _dbSlipTime->close();
  _dbPeakRate->close();

  // Allocate slip field
  _slip = new real_section_type(faultMesh->comm(), faultMesh->debug());
  _slip->setFiberDimension(vertices, spaceDim);
  faultMesh->allocate(_slip);
  assert(!_slip.isNull());
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
const ALE::Obj<pylith::real_section_type>&
pylith::faults::BruneSlipFn::slip(const double t,
				  const ALE::Obj<Mesh>& faultMesh)
{ // slip
  assert(!_parameters.isNull());
  assert(!_slip.isNull());
  assert(!faultMesh.isNull());

  const int spaceDim = _spaceDim;
  const int indexFinalSlip = 0;
  const int indexPeakRate = spaceDim + _BruneSlipFn::offsetPeakRate;
  const int indexSlipTime = spaceDim + _BruneSlipFn::offsetSlipTime;

  double_array slipValues(spaceDim);
  
  // Get vertices in fault mesh
  const ALE::Obj<Mesh::label_sequence>& vertices = faultMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();
  const int numVertices = vertices->size();

  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter) {
    const real_section_type::value_type* paramsVertex = 
      _parameters->restrictPoint(*v_iter);
    assert(0 != paramsVertex);

    const double* finalSlip = &paramsVertex[indexFinalSlip];
    const double peakRate = paramsVertex[indexPeakRate];
    const double slipTime = paramsVertex[indexSlipTime];
    
    double finalSlipMag = 0.0;
    for (int i=0; i < spaceDim; ++i)
      finalSlipMag += finalSlip[i]*finalSlip[i];
    finalSlipMag = sqrt(finalSlipMag);

    const double slip = _slipFn(t-slipTime, finalSlipMag, peakRate);
    const double scale = finalSlipMag > 0.0 ? slip / finalSlipMag : 0.0;
    for (int i=0; i < spaceDim; ++i)
      slipValues[i] = scale * finalSlip[i];

    // Update field
    _slip->updatePoint(*v_iter, &slipValues[0]);
  } // for

  PetscLogFlopsNoCheck(numVertices * (2+8 + 3*spaceDim));

  return _slip;
} // slip

// ----------------------------------------------------------------------
// Get increment of slip on fault surface between time t0 and t1.
const ALE::Obj<pylith::real_section_type>&
pylith::faults::BruneSlipFn::slipIncr(const double t0,
				      const double t1,
				      const ALE::Obj<Mesh>& faultMesh)
{ // slipIncr
  assert(!_parameters.isNull());
  assert(!_slip.isNull());
  assert(!faultMesh.isNull());

  const int spaceDim = _spaceDim;
  const int indexFinalSlip = 0;
  const int indexPeakRate = spaceDim + _BruneSlipFn::offsetPeakRate;
  const int indexSlipTime = spaceDim + _BruneSlipFn::offsetSlipTime;

  double_array slipValues(spaceDim);
  
  // Get vertices in fault mesh
  const ALE::Obj<Mesh::label_sequence>& vertices = faultMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  int count = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++count) {
    const real_section_type::value_type* paramsVertex = 
      _parameters->restrictPoint(*v_iter);
    assert(0 != paramsVertex);

    const double* finalSlip = &paramsVertex[indexFinalSlip];
    const double peakRate = paramsVertex[indexPeakRate];
    const double slipTime = paramsVertex[indexSlipTime];
    
    double finalSlipMag = 0.0;
    for (int i=0; i < spaceDim; ++i)
      finalSlipMag += finalSlip[i]*finalSlip[i];
    finalSlipMag = sqrt(finalSlipMag);

    const double slip0 = _slipFn(t0-slipTime, finalSlipMag, peakRate);
    const double slip1 = _slipFn(t1-slipTime, finalSlipMag, peakRate);
    const double scale = finalSlipMag > 0.0 ? 
      (slip1 - slip0) / finalSlipMag : 0.0;
    for (int i=0; i < spaceDim; ++i)
      slipValues[i] = scale * finalSlip[i];

    // Update field
    _slip->updatePoint(*v_iter, &slipValues[0]);
  } // for

  PetscLogFlopsNoCheck(count * (3+2*8 + 3*spaceDim));

  return _slip;
} // slipIncr

// ----------------------------------------------------------------------
// Get final slip.
ALE::Obj<pylith::real_section_type>
pylith::faults::BruneSlipFn::finalSlip(void)
{ // finalSlip
  return _parameters->getFibration(0);
} // finalSlip

// ----------------------------------------------------------------------
// Get time when slip begins at each point.
ALE::Obj<pylith::real_section_type>
pylith::faults::BruneSlipFn::slipTime(void)
{ // slipTime
  return _parameters->getFibration(2);
} // slipTime


// End of file 
