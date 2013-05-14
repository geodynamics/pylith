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

#include "TestTractPerturbation.hh" // Implementation of class methods

#include "pylith/faults/TractPerturbation.hh" // USES TractPerturbation

#include "TestFaultMesh.hh" // USES createFaultMesh()

#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestTractPerturbation );

// ----------------------------------------------------------------------
namespace pylith {
  namespace faults {
    namespace _TestTractPerturbation {
      const PylithScalar lengthScale = 1.0e+3;
      const PylithScalar pressureScale = 2.25e+10;
      const PylithScalar timeScale = 0.5;
      const PylithScalar velocityScale = lengthScale / timeScale;
      const PylithScalar densityScale = pressureScale / (velocityScale*velocityScale);
    } // namespace _TestTractPerturbation
  } // faults
} // pylith


// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestTractPerturbation::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  TractPerturbation eqsrc;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test label().
void
pylith::faults::TestTractPerturbation::testLabel(void)
{ // testLabel
  PYLITH_METHOD_BEGIN;

  const std::string& label = "nucleation";

  TractPerturbation tract;
  tract.label(label.c_str());
  CPPUNIT_ASSERT_EQUAL(label, tract._label);

  PYLITH_METHOD_END;
} // testLabel

// ----------------------------------------------------------------------
// Test hasParameter().
void
pylith::faults::TestTractPerturbation::testHasParameter(void)
{ // testHasParameter
  PYLITH_METHOD_BEGIN;

  spatialdata::spatialdb::SimpleDB db;
  
  TractPerturbation tract;

  // no values
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_initial_value"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_rate_of_change"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_change_in_value"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_rate_start_time"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_change_start_time"));

  // initial value
  tract.dbInitial(&db);
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_initial_value"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_rate_of_change"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_change_in_value"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_rate_start_time"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_change_start_time"));

  // change value
  tract.dbChange(&db);
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_initial_value"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_rate_of_change"));
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_change_in_value"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_rate_start_time"));
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_change_start_time"));

  // rate value, remove change
  tract.dbRate(&db);
  tract.dbChange(0);
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_initial_value"));
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_rate_of_change"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_change_in_value"));
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_rate_start_time"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_change_start_time"));

  PYLITH_METHOD_END;
} // testHasParameter

// ----------------------------------------------------------------------
// Test initialize() using 2-D mesh.
void
pylith::faults::TestTractPerturbation::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  topology::Mesh faultMesh;
  TractPerturbation tract;
  _initialize(&mesh, &faultMesh, &tract);
  
  // Rely on testTraction() for verification of results.

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test calculate() using 2-D mesh().
void
pylith::faults::TestTractPerturbation::testCalculate(void)
{ // testCalculate
  PYLITH_METHOD_BEGIN;

  const PylithScalar tractionE[4] = { 
    -1.0*(-2.0+1.0), -1.0*(1.0-0.5), // initial + change
    -1.0*(-2.1), -1.0*(1.1), // initial
  };

  topology::Mesh mesh;
  topology::Mesh faultMesh;
  TractPerturbation tract;
  _initialize(&mesh, &faultMesh, &tract);
  
  const PylithScalar t = 2.134 / _TestTractPerturbation::timeScale;
  tract.calculate(t);

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();CPPUNIT_ASSERT(cs);
  const int spaceDim = cs->spaceDim();

  PetscDM faultDMMesh = faultMesh.dmMesh();CPPUNIT_ASSERT(faultDMMesh);
  topology::Stratum verticesStratum(faultDMMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  CPPUNIT_ASSERT(tract._parameters);
  topology::VecVisitorMesh valueVisitor(tract._parameters->get("value"));
  const PetscScalar* valueArray = valueVisitor.localArray();CPPUNIT_ASSERT(valueArray);

  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt v = vStart, iPoint = 0; v < vEnd; ++v, ++iPoint) {
    const PetscInt voff = valueVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, valueVisitor.sectionDof(v));

    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar valueE = tractionE[iPoint*spaceDim+d];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, valueArray[voff+d]*_TestTractPerturbation::pressureScale, tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testCalculate

// ----------------------------------------------------------------------
// Test parameterFields() using 2-D mesh.
void
pylith::faults::TestTractPerturbation::testParameterFields(void)
{ // testParameterFields
  PYLITH_METHOD_BEGIN;

  const int spaceDim  = 2;
  const int fiberDimE = 7;
  const PylithScalar parametersE[2*fiberDimE] = {
    0.0, 0.0,   2.0, -1.0,   -1.0, 0.5, 1.5,
    0.0, 0.0,   2.1, -1.1,   -0.8, 0.7, 2.5,
  };
  CPPUNIT_ASSERT_EQUAL(fiberDimE, spaceDim+spaceDim+spaceDim+1);

  topology::Mesh mesh;
  topology::Mesh faultMesh;
  TractPerturbation tract;
  _initialize(&mesh, &faultMesh, &tract);
  
  const topology::Fields* parameters = tract.parameterFields();CPPUNIT_ASSERT(parameters);

  PetscDM faultDMMesh = faultMesh.dmMesh();CPPUNIT_ASSERT(faultDMMesh);
  topology::Stratum verticesStratum(faultDMMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  topology::VecVisitorMesh valueVisitor(parameters->get("value"));
  const PetscScalar* valueArray = valueVisitor.localArray();CPPUNIT_ASSERT(valueArray);

  topology::VecVisitorMesh initialVisitor(parameters->get("initial"));
  const PetscScalar* initialArray = initialVisitor.localArray();CPPUNIT_ASSERT(initialArray);

  topology::VecVisitorMesh changeVisitor(parameters->get("change"));
  const PetscScalar* changeArray = changeVisitor.localArray();CPPUNIT_ASSERT(changeArray);

  topology::VecVisitorMesh changeTimeVisitor(parameters->get("change time"));
  const PetscScalar* changeTimeArray = changeTimeVisitor.localArray();CPPUNIT_ASSERT(changeTimeArray);


  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt v = vStart, iPoint = 0; v < vEnd; ++v, ++iPoint) {
    PetscInt iE = 0;

    PetscInt off = valueVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, valueVisitor.sectionDof(v));
    for (int d=0; d < spaceDim; ++d, ++iE) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[iPoint*fiberDimE+iE], valueArray[off+d]*_TestTractPerturbation::pressureScale, tolerance);
    } // for

    off = initialVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, initialVisitor.sectionDof(v));
    for (int d=0; d < spaceDim; ++d, ++iE) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[iPoint*fiberDimE+iE], initialArray[off+d]*_TestTractPerturbation::pressureScale, tolerance);
    } // for

    off = changeVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, changeVisitor.sectionDof(v));
    for (int d=0; d < spaceDim; ++d, ++iE) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[iPoint*fiberDimE+iE], changeArray[off+d]*_TestTractPerturbation::pressureScale, tolerance);
    } // for

    off = changeTimeVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(1, changeTimeVisitor.sectionDof(v));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[iPoint*fiberDimE+iE], changeTimeArray[off]*_TestTractPerturbation::timeScale, tolerance);
    ++iE;

  } // for

  PYLITH_METHOD_END;
} // testParameterFields

// ----------------------------------------------------------------------
// Test vertexField() using 2-D mesh.
void
pylith::faults::TestTractPerturbation::testVertexField(void)
{ // testVertexField
  PYLITH_METHOD_BEGIN;

  const int fiberDimE = 1;
  const PylithScalar fieldE[2*fiberDimE] = {
    1.5,
    2.5,
  };
  const char* label = "traction_change_start_time";

  topology::Mesh mesh;
  topology::Mesh faultMesh;
  TractPerturbation tract;
  _initialize(&mesh, &faultMesh, &tract);

  const topology::Field& field = tract.vertexField(label);

  PetscDM faultDMMesh = faultMesh.dmMesh();CPPUNIT_ASSERT(faultDMMesh);
  topology::Stratum verticesStratum(faultDMMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  topology::VecVisitorMesh fieldVisitor(field);
  const PetscScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);

  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt v = vStart, iPoint = 0; v < vEnd; ++v, ++iPoint) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, fieldVisitor.sectionDof(v));

    for(PetscInt d = 0; d < fiberDimE; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(fieldE[iPoint*fiberDimE+d], fieldArray[off+d]*_TestTractPerturbation::timeScale, tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testVertexField

// ----------------------------------------------------------------------
// Initialize TractPerturbation.
void
pylith::faults::TestTractPerturbation::_initialize(topology::Mesh* mesh,
						   topology::Mesh* faultMesh,
						   TractPerturbation* tract)
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(faultMesh);
  CPPUNIT_ASSERT(tract);
  PetscErrorCode err;

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* initialFilename = "data/tri3_initialtract.spatialdb";
  const char* changeFilename = "data/tri3_changetract.spatialdb";
  const PylithScalar orientationVertex[4] = {
      0.0, -1.0,  -1.0, 0.0,
  };

  meshio::MeshIOAscii meshIO;
  meshIO.filename(meshFilename);
  meshIO.debug(false);
  meshIO.interpolate(false);
  meshIO.read(mesh);

  // Set coordinate system
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  const int spaceDim = cs.spaceDim();

  // Set scales
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_TestTractPerturbation::lengthScale);
  normalizer.pressureScale(_TestTractPerturbation::pressureScale);
  normalizer.densityScale(_TestTractPerturbation::densityScale);
  normalizer.timeScale(_TestTractPerturbation::timeScale);
  topology::MeshOps::nondimensionalize(mesh, normalizer);

  // Create fault mesh
  TestFaultMesh::createFaultMesh(faultMesh, mesh, faultLabel, faultId);

  // Setup databases
  spatialdata::spatialdb::SimpleDB dbInitial("initial traction");
  spatialdata::spatialdb::SimpleIOAscii ioInitial;
  ioInitial.filename(initialFilename);
  dbInitial.ioHandler(&ioInitial);
  
  spatialdata::spatialdb::SimpleDB dbChange("traction change");
  spatialdata::spatialdb::SimpleIOAscii ioChange;
  ioChange.filename(changeFilename);
  dbChange.ioHandler(&ioChange);

  // Setup fault orientation
  topology::Field faultOrientation(*faultMesh);
  faultOrientation.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim*spaceDim);
  faultOrientation.allocate();

  PetscDM faultDMMesh = faultMesh->dmMesh();CPPUNIT_ASSERT(faultDMMesh);
  topology::Stratum verticesStratum(faultDMMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  topology::VecVisitorMesh orientationVisitor(faultOrientation);
  PetscScalar* orientationArray = orientationVisitor.localArray();CPPUNIT_ASSERT(orientationArray);

  const PetscInt orientationSize = spaceDim*spaceDim;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt off = orientationVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(orientationSize, orientationVisitor.sectionDof(v));
    
    for(PetscInt d = 0; d < orientationSize; ++d) {
      orientationArray[off+d] = orientationVertex[d];
    }
  } // for
  
  // setup TractPerturbation
  tract->dbInitial(&dbInitial);
  tract->dbChange(&dbChange);
  
  tract->initialize(*faultMesh, faultOrientation, normalizer);

  PYLITH_METHOD_END;
} // _initialize


// End of file 
