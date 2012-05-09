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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestTractPerturbation.hh" // Implementation of class methods

#include "pylith/faults/TractPerturbation.hh" // USES TractPerturbation

#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldsNew.hh" // USES FieldsNew
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestTractPerturbation );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::RealUniformSection RealUniformSection;

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestTractPerturbation::testConstructor(void)
{ // testConstructor
  TractPerturbation eqsrc;
} // testConstructor

// ----------------------------------------------------------------------
// Test label().
void
pylith::faults::TestTractPerturbation::testLabel(void)
{ // testLabel
  const std::string& label = "nucleation";

  TractPerturbation tract;
  tract.label(label.c_str());
  CPPUNIT_ASSERT_EQUAL(label, tract._label);
} // testLabel

// ----------------------------------------------------------------------
// Test hasParameter().
void
pylith::faults::TestTractPerturbation::testHasParameter(void)
{ // testHasParameter
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
} // testHasParameter

// ----------------------------------------------------------------------
// Test initialize() using 2-D mesh.
void
pylith::faults::TestTractPerturbation::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  TractPerturbation tract;
  _initialize(&mesh, &faultMesh, &tract);
  
  // Rely on testTraction() for verification of results.
} // testInitialize

// ----------------------------------------------------------------------
// Test calculate() using 2-D mesh().
void
pylith::faults::TestTractPerturbation::testCalculate(void)
{ // testCalculate
  const PylithScalar tractionE[4] = { 
    -1.0*(-2.0+1.0), -1.0*(1.0-0.5), // initial + change
    -1.0*(-2.1), -1.0*(1.1), // initial
  };

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  TractPerturbation tract;
  _initialize(&mesh, &faultMesh, &tract);
  
  const PylithScalar t = 2.134;
  tract.calculate(t);


  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  CPPUNIT_ASSERT(cs);
  const int spaceDim = cs->spaceDim();

  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh.sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = faultSieveMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

  //traction.view("TRACTION"); // DEBUGGING

  const PylithScalar tolerance = 1.0e-06;
  int iPoint = 0;

  CPPUNIT_ASSERT(tract._parameters);
  const ALE::Obj<RealUniformSection>& paramsSection = tract._parameters->section();
  CPPUNIT_ASSERT(!paramsSection.isNull());
  const int fiberDim = tract._parameters->fiberDim();
  const int valueIndex = tract._parameters->sectionIndex("value");
  const int valueFiberDim = tract._parameters->sectionFiberDim("value");
  CPPUNIT_ASSERT_EQUAL(spaceDim, valueFiberDim);
  CPPUNIT_ASSERT(valueIndex + valueFiberDim < fiberDim);
  
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin(); v_iter != verticesEnd; ++v_iter, ++iPoint) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, paramsSection->getFiberDimension(*v_iter));
    const PylithScalar* vals = paramsSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(vals);
    
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      const PylithScalar valueE = tractionE[iPoint*spaceDim+iDim];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, vals[valueIndex+iDim], tolerance);
    } // for
  } // for
} // testCalculate

// ----------------------------------------------------------------------
// Test parameterFields() using 2-D mesh.
void
pylith::faults::TestTractPerturbation::testParameterFields(void)
{ // testParameterFields
  const int fiberDimE = 7;
  const PylithScalar parametersE[2*fiberDimE] = {
    0.0, 0.0,   2.0, -1.0,   -1.0, 0.5, 1.5,
    0.0, 0.0,   2.1, -1.1,   -0.8, 0.7, 2.5,
  };

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  TractPerturbation tract;
  _initialize(&mesh, &faultMesh, &tract);
  
  const topology::FieldsNew<topology::SubMesh>* parameters = tract.parameterFields();
  const ALE::Obj<RealUniformSection>& parametersSection = parameters->section();
  CPPUNIT_ASSERT(!parametersSection.isNull());

  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh.sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = faultSieveMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

  const PylithScalar tolerance = 1.0e-06;
  int iPoint = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin(); v_iter != verticesEnd; ++v_iter, ++iPoint) {
    const int fiberDim = parametersSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
    const PylithScalar* vals = parametersSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(vals);

    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const PylithScalar valueE = parametersE[iPoint*fiberDim+iDim];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, vals[iDim], tolerance);
    } // for
  } // for
} // testParameterFields

// ----------------------------------------------------------------------
// Test vertexField() using 2-D mesh.
void
pylith::faults::TestTractPerturbation::testVertexField(void)
{ // testVertexField
  const int fiberDimE = 1;
  const PylithScalar fieldE[2*fiberDimE] = {
    1.5,
    2.5,
  };
  const char* label = "traction_change_start_time";

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  TractPerturbation tract;
  _initialize(&mesh, &faultMesh, &tract);

  const topology::Field<topology::SubMesh>& field = tract.vertexField(label);
  const ALE::Obj<RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh.sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = faultSieveMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

  const PylithScalar tolerance = 1.0e-06;
  int iPoint = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin(); v_iter != verticesEnd; ++v_iter, ++iPoint) {
    const int fiberDim = section->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
    const PylithScalar* vals = section->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(vals);

    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const PylithScalar valueE = fieldE[iPoint*fiberDim+iDim];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, vals[iDim], tolerance);
    } // for
  } // for
} // testVertexField

// ----------------------------------------------------------------------
// Initialize TractPerturbation.
void
pylith::faults::TestTractPerturbation::_initialize(topology::Mesh* mesh,
						   topology::SubMesh* faultMesh,
						   TractPerturbation* tract)
{ // _initialize
  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(faultMesh);
  CPPUNIT_ASSERT(tract);

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

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  const int spaceDim = cs.spaceDim();

  // Create fault mesh
  int firstFaultVertex = 0;
  int firstLagrangeVertex = mesh->sieveMesh()->getIntSection(faultLabel)->size();
  int firstFaultCell = mesh->sieveMesh()->getIntSection(faultLabel)->size();
  const bool useLagrangeConstraints = true;
  if (useLagrangeConstraints) {
    firstFaultCell += mesh->sieveMesh()->getIntSection(faultLabel)->size();
  } // if
  ALE::Obj<SieveFlexMesh> faultBoundary = 0;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CohesiveTopology::createFault(faultMesh, faultBoundary,
                                *mesh, sieveMesh->getIntSection(faultLabel));
  CohesiveTopology::create(mesh, *faultMesh, faultBoundary, 
                           sieveMesh->getIntSection(faultLabel),
                           faultId,
                           firstFaultVertex, firstLagrangeVertex, firstFaultCell,
                           useLagrangeConstraints);
  // Need to copy coordinates from mesh to fault mesh since we are
  // using create() instead of createParallel().
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh->sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  faultSieveMesh->setRealSection("coordinates", 
				 sieveMesh->getRealSection("coordinates"));

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
  const ALE::Obj<SieveMesh::label_sequence>& vertices = faultSieveMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  topology::Field<topology::SubMesh> faultOrientation(*faultMesh);
  faultOrientation.newSection(vertices, spaceDim*spaceDim);
  faultOrientation.allocate();
  const ALE::Obj<RealSection>& orientationSection = faultOrientation.section();
  CPPUNIT_ASSERT(!orientationSection.isNull());
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin(); v_iter != verticesEnd; ++v_iter) {
    CPPUNIT_ASSERT_EQUAL(spaceDim*spaceDim, orientationSection->getFiberDimension(*v_iter));
    orientationSection->updatePoint(*v_iter, orientationVertex);
  } // for
  
  spatialdata::units::Nondimensional normalizer;

  // setup TractPerturbation
  tract->dbInitial(&dbInitial);
  tract->dbChange(&dbChange);
  
  tract->initialize(*faultMesh, faultOrientation, normalizer);
} // _initialize


// End of file 
