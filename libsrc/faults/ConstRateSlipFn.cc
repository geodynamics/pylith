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

#include "ConstRateSlipFn.hh" // implementation of object methods

#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#include "pylith/utils/array.hh" // USES double_array

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

namespace pylith {
  namespace faults {
    namespace _ConstRateSlipFn {
      const int offsetSlipTime = 0;
    } // _ConstRateSlipFn
  } // faults
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::ConstRateSlipFn::ConstRateSlipFn(void) :
  _dbSlipRate(0),
  _dbSlipTime(0),
  _spaceDim(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::ConstRateSlipFn::~ConstRateSlipFn(void)
{ // destructor
  _dbSlipRate = 0;
  _dbSlipTime = 0;
} // destructor

// ----------------------------------------------------------------------
// Initialize slip time function.
void
pylith::faults::ConstRateSlipFn::initialize(
			   const ALE::Obj<Mesh>& faultMesh,
			   const spatialdata::geocoords::CoordSys* cs,
			   const spatialdata::units::Nondimensional& normalizer,
			   const double originTime)
{ // initialize
  assert(!faultMesh.isNull());
  assert(0 != cs);
  assert(0 != _dbSlipRate);
  assert(0 != _dbSlipTime);

  _spaceDim = cs->spaceDim();
  const int spaceDim = _spaceDim;
  const int indexSlipRate = 0;
  const int indexSlipTime = spaceDim + _ConstRateSlipFn::offsetSlipTime;

  // Get vertices in fault mesh
  const ALE::Obj<Mesh::label_sequence>& vertices = faultMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  const int fiberDim = spaceDim + 1;
  _parameters = new real_section_type(faultMesh->comm(), faultMesh->debug());
  _parameters->addSpace(); // slip rate
  _parameters->addSpace(); // slip time
  assert(2 == _parameters->getNumSpaces());
  _parameters->setChart(real_section_type::chart_type(*std::min_element(vertices->begin(), vertices->end()), *std::max_element(vertices->begin(), vertices->end())+1));
  _parameters->setFiberDimension(vertices, fiberDim);
  _parameters->setFiberDimension(vertices, spaceDim, 0); // slip rate
  _parameters->setFiberDimension(vertices, 1, 1); // slip time
  faultMesh->allocate(_parameters);
  assert(!_parameters.isNull());

  // Open databases and set query values
  _dbSlipRate->open();
  switch (spaceDim)
    { // switch
    case 1 : {
      const char* slipValues[] = {"fault-opening"};
      _dbSlipRate->queryVals(slipValues, 1);
      break;
    } // case 1
    case 2 : {
      const char* slipValues[] = {"left-lateral-slip", "fault-opening"};
      _dbSlipRate->queryVals(slipValues, 2);
      break;
    } // case 2
    case 3 : {
      const char* slipValues[] = {"left-lateral-slip", "reverse-slip", 
				  "fault-opening"};
      _dbSlipRate->queryVals(slipValues, 3);
      break;
    } // case 3
    default :
      assert(0);
    } // switch

  _dbSlipTime->open();
  const char* slipTimeValues[] = {"slip-time"};
  _dbSlipTime->queryVals(slipTimeValues, 1);

  // Get coordinates of vertices
  const ALE::Obj<real_section_type>& coordinates = 
    faultMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  const double lengthScale = normalizer.lengthScale();
  const double timeScale = normalizer.timeScale();
  const double velocityScale =
    normalizer.lengthScale() / normalizer.timeScale();

  double_array paramsVertex(fiberDim);
  double_array vCoordsGlobal(spaceDim);
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter) {
    coordinates->restrictPoint(*v_iter, 
			       &vCoordsGlobal[0], vCoordsGlobal.size());
    normalizer.dimensionalize(&vCoordsGlobal[0], vCoordsGlobal.size(),
			      lengthScale);
    
    int err = _dbSlipRate->query(&paramsVertex[indexSlipRate], spaceDim, 
				 &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find slip rate at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbSlipRate->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&paramsVertex[indexSlipRate], spaceDim,
				 velocityScale);

    err = _dbSlipTime->query(&paramsVertex[indexSlipTime], 1, 
			     &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find slip initiation time at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbSlipTime->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&paramsVertex[indexSlipTime], 1,
				 timeScale);
    // add origin time to rupture time
    paramsVertex[indexSlipTime] += originTime;

    _parameters->updatePoint(*v_iter, &paramsVertex[0]);
  } // for

  // Close databases
  _dbSlipRate->close();
  _dbSlipTime->close();
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
void
pylith::faults::ConstRateSlipFn::slip(const ALE::Obj<pylith::real_section_type>& slipField,
				      const double t,
				      const ALE::Obj<Mesh>& faultMesh)
{ // slip
  assert(!_parameters.isNull());
  assert(!slipField.isNull());
  assert(!faultMesh.isNull());

  const int spaceDim = _spaceDim;
  const int indexSlipRate = 0;
  const int indexSlipTime = spaceDim + _ConstRateSlipFn::offsetSlipTime;

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

    const double* slipRate = &paramsVertex[indexSlipRate];
    const double slipTime = paramsVertex[indexSlipTime];
    slipValues = 0.0;
    
    const double relTime = t - slipTime;
    if (relTime > 0)
      for (int i=0; i < spaceDim; ++i)
	slipValues[i] = slipRate[i] * relTime;
    
    // Update field
    slipField->updateAddPoint(*v_iter, &slipValues[0]);
  } // for

  PetscLogFlops(numVertices * (1+1 + 4*spaceDim));
} // slip

// ----------------------------------------------------------------------
// Get increment of slip on fault surface between time t0 and t1.
void
pylith::faults::ConstRateSlipFn::slipIncr(const ALE::Obj<pylith::real_section_type>& slipField,
					  const double t0,
					  const double t1,
					  const ALE::Obj<Mesh>& faultMesh)
{ // slipIncr
  assert(!_parameters.isNull());
  assert(!slipField.isNull());
  assert(!faultMesh.isNull());

  const int spaceDim = _spaceDim;
  const int indexSlipRate = 0;
  const int indexSlipTime = spaceDim + _ConstRateSlipFn::offsetSlipTime;

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

    const double* slipRate = &paramsVertex[indexSlipRate];
    const double slipTime = paramsVertex[indexSlipTime];
    slipValues = 0.0;

    const double relTime0 = t0 - slipTime;
    const double relTime1 = t1 - slipTime;
    double elapsedTime = 0.0;
    if (relTime0 > 0)
      elapsedTime = t1 - t0;
    else if (relTime1 > 0)
      elapsedTime = t1 - slipTime;
    for (int i=0; i < spaceDim; ++i)
      slipValues[i] = slipRate[i] * elapsedTime;
    
    // Update field
    slipField->updateAddPoint(*v_iter, &slipValues[0]);
  } // for

  PetscLogFlops(numVertices * (2 + spaceDim));
} // slipIncr

// ----------------------------------------------------------------------
// Get final slip (slip rate in this case).
ALE::Obj<pylith::real_section_type>
pylith::faults::ConstRateSlipFn::finalSlip(void)
{ // finalSlip
  // This is actually slip rate.
  return _parameters->getFibration(0);
} // finalSlip

// ----------------------------------------------------------------------
// Get time when slip begins at each point.
ALE::Obj<pylith::real_section_type>
pylith::faults::ConstRateSlipFn::slipTime(void)
{ // slipTime
  return _parameters->getFibration(1);
} // slipTime


// End of file 
