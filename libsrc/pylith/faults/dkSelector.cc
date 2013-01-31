// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
// Romain Jolivet, Caltech
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "dkSelector.hh" // implementation of object methods

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
pylith::faults::dkSelector::dkSelector(void) :
  _dbdksel(0),
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::dkSelector::~dkSelector(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::dkSelector::deallocate(void)
{ // deallocate
  _dbdksel = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize dkselector
void
pylith::faults::dkSelector::initialize(
			    const topology::SubMesh& faultMesh,
			    const spatialdata::units::Nondimensional& normalizer)
{ // initialize
  assert(0 != _dbdksel);

  // Get the spatial coordinate and the dimension of the problem
  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  // Get the normalizing values
  const PylithScalar lengthScale = normalizer.lengthScale();

  // Get vertices in fault mesh
  const ALE::Obj<SieveMesh>& sieveMesh = faultMesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const label_sequence::iterator verticesBegin = vertices->begin();
  const label_sequence::iterator verticesEnd = vertices->end();

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Fault");

  // Create a parameter to go fetch into the spatial data (DKsel is a par array between 0 and 1; >0.5 is kinematic, <=0.5 is dynamic)
  delete _parameters; _parameters = new topology::Fields<topology::Field<topology::SubMesh> >(faultMesh);
  assert(0 != _parameters);
  _parameters->add("Dynamic Kinematic Selector","dynamic_kinematic_selector");
  topology::Field<topology::SubMesh>& DKSel = _parameters->get("Dynamic Kinematic Selector");
  DKSel.newSection(vertices, spaceDim);
  DKsel.allocate();
  DKsel.vectorFieldType(topology::FieldBase::VECTOR);
  const ALE::Obj<RealSection>& DKselSection = DKsel.section();
  assert(!DKselSection.isNull());  

  logger.stagePop();

  // Open databases and set query values
  _dbdksel->open();
  const char* dkselValues[] = {"dynamic-kinematic"};
  _dbdksel->queryVals(dkselValues, 1);

  // Get coordinates of vertices
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  _dkselVertex.resize(spaceDim);
  scalar_array vCoordsGlobal(spaceDim);
  for (label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    coordinates->restrictPoint(*v_iter, 
			       &vCoordsGlobal[0], vCoordsGlobal.size());
    normalizer.dimensionalize(&vCoordsGlobal[0], vCoordsGlobal.size(),
			      lengthScale);
    
    int err = _dbdksel->query(&_dkselVertex[0], _dkSelVertex.size(), 
				 &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find DK Selector at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbdksel->label() << ".";
      throw std::runtime_error(msg.str());
    } // if

    DKselSection->updatePoint(*v_iter, &_dkselVertex[0]);
  } // for

  // Close databases
  _dbdksel->close();
} // initialize

// ----------------------------------------------------------------------
// Get dynamic kinematic selector field on fault surface
void
pylith::faults::dkSelector::dk(topology::Field<topology::SubMesh>* dk)
{ // dk
  assert(0 != dk);
  assert(0 != _parameters);

  // Get vertices in fault mesh
  const ALE::Obj<SieveMesh>& sieveMesh = dk->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const label_sequence::iterator verticesBegin = vertices->begin();
  const label_sequence::iterator verticesEnd = vertices->end();

  // Get sections
  const topology::Field<topology::SubMesh>& DKSel = 
    _parameters->get("Dynamic Kinematic Selector");
  const ALE::Obj<RealSection>& DKselSection = DKsel.section();
  assert(!DKselSection.isNull());

  for (label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    DKselSection->restrictPoint(*v_iter, &_dkselVertex[0],
				   _dkselVertex.size());
    // if dkselVertex is under 0.5, the vertex is kinematically controled
    if ( _dkselVertex < 0.5 ){
	_dkselVertex = 0.0;
    else {
	_dkselVertex = 1.0;
    }
    DKselSection->updateAddPoint(*v_iter, &_dkselVertex[0]);
  } // for

  PetscLogFlops(vertices->size());
} // slip

// End of file 
