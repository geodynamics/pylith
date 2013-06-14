// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "StaticPerturbation.hh" // implementation of object methods

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveMesh;
typedef pylith::topology::SubMesh::SieveMesh::label_sequence label_sequence;
typedef pylith::topology::SubMesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::StaticPerturbation::StaticPerturbation(void) :
  _dbAmplitude(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::StaticPerturbation::~StaticPerturbation(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::StaticPerturbation::deallocate(void)
{ // deallocate
  TractPerturbation::deallocate();

  _dbAmplitude = 0; // :TODO: Use shared pointer
} // deallocate
  
// ----------------------------------------------------------------------
// Set spatial database for amplitude of traction.
void
pylith::faults::StaticPerturbation::dbAmplitude(spatialdata::spatialdb::SpatialDB* const db) {
  _dbAmplitude = db;
} // dbAmplitude

// ----------------------------------------------------------------------
// Initialize traction perturbation function.
void
pylith::faults::StaticPerturbation::initialize(
			    const topology::SubMesh& faultMesh,
			    const spatialdata::units::Nondimensional& normalizer)
{ // initialize
  assert(_dbAmplitude);

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  const PylithScalar lengthScale = normalizer.lengthScale();
  const PylithScalar timeScale = normalizer.timeScale();
  const PylithScalar pressureScale = normalizer.pressureScale();

  // Memory logging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Fault");

  // Get vertices in fault mesh
  const ALE::Obj<SieveMesh>& sieveMesh = faultMesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const label_sequence::iterator verticesBegin = vertices->begin();
  const label_sequence::iterator verticesEnd = vertices->end();

  delete _parameters; _parameters = new topology::Fields<topology::Field<topology::SubMesh> >(faultMesh);
  assert(_parameters);
  _parameters->add("amplitude", "amplitude");

  topology::Field<topology::SubMesh>& amplitude = _parameters->get("amplitude");
  amplitude.newSection(vertices, spaceDim);
  amplitude.allocate();
  amplitude.scale(pressureScale);
  amplitude.vectorFieldType(topology::FieldBase::VECTOR);
  const ALE::Obj<RealSection>& amplitudeSection = amplitude.section();
  assert(!amplitudeSection.isNull());  

  logger.stagePop();

  // Open databases and set query values
  _dbAmplitude->open();
  switch (spaceDim)
    { // switch
    case 1 : {
      const char* tractionValues[1] = {"traction-normal"};
      _dbAmplitude->queryVals(tractionValues, 1);
      break;
    } // case 1
    case 2 : {
      const char* tractionValues[2] = {"traction-shear", "traction-normal"};
      _dbAmplitude->queryVals(tractionValues, 2);
      break;
    } // case 2
    case 3 : {
      const char* tractionValues[3] = {"traction-shear-leftlateral", 
				   "traction-shear-updip",
				   "traction-normal"};
      _dbAmplitude->queryVals(tractionValues, 3);
      break;
    } // case 3
    default :
      std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad spatial dimension in StaticPerturbation.");
    } // switch

  // Get coordinates of vertices
  const ALE::Obj<RealSection>& coordinates = sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  scalar_array tractionVertex(spaceDim);
  scalar_array vCoordsGlobal(spaceDim);  
  for (label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    coordinates->restrictPoint(*v_iter, &vCoordsGlobal[0], vCoordsGlobal.size());
    normalizer.dimensionalize(&vCoordsGlobal[0], vCoordsGlobal.size(), lengthScale);
        
    int err = _dbAmplitude->query(&tractionVertex[0], tractionVertex.size(), 
				  &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find traction at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbAmplitude->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&tractionVertex[0], tractionVertex.size(), pressureScale);

    assert(spaceDim == amplitudeSection->getFiberDimension(*v_iter));
    amplitudeSection->updatePoint(*v_iter, &tractionVertex[0]);
  } // for

  // Close databases
  _dbAmplitude->close();
} // initialize

// ----------------------------------------------------------------------
// Get traction on fault surface at time t.
void
pylith::faults::StaticPerturbation::traction(topology::Field<topology::SubMesh>* tractionField,
				 const PylithScalar t)
{ // traction
  assert(tractionField);
  assert(_parameters);

  const spatialdata::geocoords::CoordSys* cs = tractionField->mesh().coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  // Get vertices in fault mesh
  const ALE::Obj<SieveMesh>& sieveMesh = tractionField->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const label_sequence::iterator verticesBegin = vertices->begin();
  const label_sequence::iterator verticesEnd = vertices->end();

  // Get sections
  const ALE::Obj<RealSection>& amplitudeSection = _parameters->get("amplitude").section();
  assert(!amplitudeSection.isNull());
  const ALE::Obj<RealSection>& tractionSection = tractionField->section();
  assert(!tractionSection.isNull());

  for (label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    assert(spaceDim == amplitudeSection->getFiberDimension(*v_iter));
    const PylithScalar* amplitudeVertex = amplitudeSection->restrictPoint(*v_iter);

    // Update field
    assert(spaceDim == tractionSection->getFiberDimension(*v_iter));
    tractionSection->updateAddPoint(*v_iter, &amplitudeVertex[0]);
  } // for

} // traction

// ----------------------------------------------------------------------
// Get traction amplitude..
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::StaticPerturbation::amplitude(void)
{ // amplitude
  return _parameters->get("amplitude");
} // amplitude


// End of file 
