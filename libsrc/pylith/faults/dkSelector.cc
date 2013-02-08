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

#include "DKSelector.hh" // implementation of object methods

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
typedef pylith::topology::SubMesh::RealUniformSection RealUniformSection;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::DKSelector::DKSelector(void) :
  _dbdksel(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::DKSelector::~DKSelector(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::DKSelector::deallocate(void)
{ // deallocate
  _dbdksel = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize dkselector
void
pylith::faults::DKSelector::initialize(
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

  // Use parameter to go fetch into the spatial data (DKSel is a par array between 0 and 1; >0.5 is kinematic, <=0.5 is dynamic)
  delete _parameters;
  _parameters = new topology::Fields<topology::Field<topology::SubMesh> >(faultMesh);
  _parameters->add("Dynamic Kinematic Selector","dynamic_kinematic_selector");
  topology::Field<topology::SubMesh>& DKSel = _parameters->get("Dynamic Kinematic Selector");
  DKSel.newSection(vertices, spaceDim);
  DKSel.allocate();
  DKSel.vectorFieldType(topology::FieldBase::VECTOR);
  const ALE::Obj<RealSection>& DKSelSection = DKSel.section();
  assert(!DKSelSection.isNull());  

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
    
    int err = _dbdksel->query(&_dkselVertex[0], _dkselVertex.size(), 
				 &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find DK Selector at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbdksel->label() << ".";
      throw std::runtime_error(msg.str());
    } // if

    DKSelSection->updatePoint(*v_iter, &_dkselVertex[0]);
  } // for

  // Close databases
  _dbdksel->close();
} // initialize

// ----------------------------------------------------------------------
// Get dynamic kinematic selector field on fault surface (time will be the argument in the future)
void
pylith::faults::DKSelector::dk(topology::Field<topology::SubMesh>* const dk)
{ // dk
  assert(0 != _parameters);

  // Get fault mesh
  const ALE::Obj<SieveMesh>& sieveMesh = dk->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const label_sequence::iterator verticesBegin = vertices->begin();
  const label_sequence::iterator verticesEnd = vertices->end();

  // Build the section
  const topology::Field<topology::SubMesh>& DKSel = _parameters->get("Dynamic Kinematic Selector");
  const ALE::Obj<RealSection>& DKSelSection = DKSel.section();
  assert(!DKSelSection.isNull());
  const ALE::Obj<RealSection>& dkSection = dk->section();

  // Iterate over vertices
  for (label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
 
    // take the good vertex in DKSelSection, put it in _dkselvv
    DKSelSection->restrictPoint(*v_iter, &_dkselvv,1); 
    // if _dkselvv is under 0.5, the vertex is kinematically controled (put the time condition here)
    if ( _dkselvv< 0.5 ){
	_dkselvv = 0.0;
    } else {
	_dkselvv = 1.0;
    }
    // Put that thing in dkSection
    dkSection->updatePoint(*v_iter, &_dkselvv);

  } // for

  PetscLogFlops(vertices->size());
} // dk

// End of file 
