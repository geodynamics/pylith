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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestStepSlipFn.hh" // Implementation of class methods

#include "pylith/faults/StepSlipFn.hh" // USES StepSlipFn

#include "TestFaultMesh.hh" // USES createFaultMesh()

#include "pylith/topology/Mesh.hh" // USES Mesh
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
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestStepSlipFn );

// ----------------------------------------------------------------------
namespace pylith {
  namespace faults {
    namespace _TestStepSlipFn {
      struct DataStruct {
	const char* meshFilename;
	const char* faultLabel;
	const int faultId;
	const char* finalSlipFilename;
	const char* slipTimeFilename;
	const int* constraintPts;
	const PylithScalar* finalSlipE;
	const PylithScalar* slipTimeE;
	const int numConstraintPts;
      }; // DataStruct
    } // _TestStepSlipFn
  } // faults
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestStepSlipFn::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  StepSlipFn slipfn;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbFinalSlip().
void
pylith::faults::TestStepSlipFn::testDbFinalSlip(void)
{ // testDbFinalSlip
  PYLITH_METHOD_BEGIN;

  const char* label = "database ABC";
  StepSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbFinalSlip(&db);

  CPPUNIT_ASSERT(slipfn._dbFinalSlip);
  CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(slipfn._dbFinalSlip->label()));
  CPPUNIT_ASSERT(!slipfn._dbSlipTime);

  PYLITH_METHOD_END;
} // testDbFinalSlip

// ----------------------------------------------------------------------
// Test dbSlipTime().
void
pylith::faults::TestStepSlipFn::testDbSlipTime(void)
{ // testDbSlipTime
  PYLITH_METHOD_BEGIN;

  const char* label = "database ABCD";
  StepSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbSlipTime(&db);

  CPPUNIT_ASSERT(slipfn._dbSlipTime);
  CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(slipfn._dbSlipTime->label()));
  CPPUNIT_ASSERT(!slipfn._dbFinalSlip);

  PYLITH_METHOD_END;
} // testDbSlipTime

// ----------------------------------------------------------------------
// Test initialize() in 2-D.
void
pylith::faults::TestStepSlipFn::testInitialize2D(void)
{ // testInitialize2D
  PYLITH_METHOD_BEGIN;

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const int constraintPts[] = { 3, 4 };
  const PylithScalar finalSlipE[] = { 2.3, 0.1, 
				2.4, 0.2 };
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const int numConstraintPts = 2;

  _TestStepSlipFn::DataStruct data = {meshFilename,
				      faultLabel,
				      faultId,
				      finalSlipFilename,
				      slipTimeFilename,
				      constraintPts,
				      finalSlipE,
				      slipTimeE,
				      numConstraintPts};
  _testInitialize(data);

  PYLITH_METHOD_END;
} // testInitialize2D

// ----------------------------------------------------------------------
// Test initialize() in 3-D.
void
pylith::faults::TestStepSlipFn::testInitialize3D(void)
{ // testInitialize3D
  PYLITH_METHOD_BEGIN;

  const char* meshFilename = "data/tet4.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tet4_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tet4_sliptime.spatialdb";
  const int constraintPts[] = { 3, 4, 5 };
  const PylithScalar finalSlipE[] = { 2.3, -0.7, 0.1,
				2.4, -0.8, 0.2,
				2.5, -0.9, 0.3 };
  const PylithScalar slipTimeE[] = { 1.2, 1.3, 1.4 };
  const int numConstraintPts = 3;

  _TestStepSlipFn::DataStruct data = {meshFilename,
				      faultLabel,
				      faultId,
				      finalSlipFilename,
				      slipTimeFilename,
				      constraintPts,
				      finalSlipE,
				      slipTimeE,
				      numConstraintPts};
  _testInitialize(data);

  PYLITH_METHOD_END;
} // testInitialize3D

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestStepSlipFn::testSlip(void)
{ // testSlip
  PYLITH_METHOD_BEGIN;

  const PylithScalar slipE[] = { 2.3, 0.1, 
			   0.0, 0.0};
  const PylithScalar originTime = 5.064;

  topology::Mesh mesh;
  topology::Mesh faultMesh;
  StepSlipFn slipfn;
  _initialize(&mesh, &faultMesh, &slipfn, originTime);
  
  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();CPPUNIT_ASSERT(cs);
  const int spaceDim = cs->spaceDim();

  topology::Field slip(faultMesh);
  slip.newSection(topology::Field::VERTICES_FIELD, spaceDim);
  slip.allocate();

  const PylithScalar t = 1.234;
  slipfn.slip(&slip, originTime+t);

  PetscDM dmMesh = faultMesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  topology::VecVisitorMesh slipVisitor(slip);
  const PetscScalar* slipArray = slipVisitor.localArray();CPPUNIT_ASSERT(slipArray);

  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt v = vStart, iPoint = 0; v < vEnd; ++v, ++iPoint) {

    const PetscInt off = slipVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, slipVisitor.sectionDof(v));

    for(PetscInt d = 0; d < spaceDim; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE[iPoint*spaceDim+d], slipArray[off+d], tolerance);
  } // for

  PYLITH_METHOD_END;
} // testSlip

// ----------------------------------------------------------------------
// Initialize StepSlipFn.
void
pylith::faults::TestStepSlipFn::_initialize(topology::Mesh* mesh,
					    topology::Mesh* faultMesh,
					    StepSlipFn* slipfn,
					    const PylithScalar originTime)
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(faultMesh);
  CPPUNIT_ASSERT(slipfn);
  PetscErrorCode err;

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";

  meshio::MeshIOAscii meshIO;
  meshIO.filename(meshFilename);
  meshIO.debug(false);
  meshIO.interpolate(false);
  meshIO.read(mesh);

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  const int spaceDim = mesh->dimension();
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);

  // Create fault mesh
  TestFaultMesh::createFaultMesh(faultMesh, mesh, faultLabel, faultId);

  // Setup databases
  spatialdata::spatialdb::SimpleDB dbFinalSlip("final slip");
  spatialdata::spatialdb::SimpleIOAscii ioFinalSlip;
  ioFinalSlip.filename(finalSlipFilename);
  dbFinalSlip.ioHandler(&ioFinalSlip);
  
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
  ioSlipTime.filename(slipTimeFilename);
  dbSlipTime.ioHandler(&ioSlipTime);
  
  spatialdata::units::Nondimensional normalizer;

  // setup StepSlipFn
  slipfn->dbFinalSlip(&dbFinalSlip);
  slipfn->dbSlipTime(&dbSlipTime);
  
  slipfn->initialize(*faultMesh, normalizer, originTime);

  PYLITH_METHOD_END;
} // _initialize

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestStepSlipFn::_testInitialize(const _TestStepSlipFn::DataStruct& data)
{ // _testInitialize
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err;
  // Setup mesh
  topology::Mesh mesh;
  meshio::MeshIOAscii meshIO;
  meshIO.filename(data.meshFilename);
  meshIO.debug(false);
  meshIO.interpolate(false);
  meshIO.read(&mesh);

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  const int spaceDim = mesh.dimension();
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh.coordsys(&cs);

  // Create fault mesh
  topology::Mesh faultMesh;
  TestFaultMesh::createFaultMesh(&faultMesh, &mesh, data.faultLabel, data.faultId);

  // Setup databases
  spatialdata::spatialdb::SimpleDB dbFinalSlip("final slip");
  spatialdata::spatialdb::SimpleIOAscii ioFinalSlip;
  ioFinalSlip.filename(data.finalSlipFilename);
  dbFinalSlip.ioHandler(&ioFinalSlip);
  
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
  ioSlipTime.filename(data.slipTimeFilename);
  dbSlipTime.ioHandler(&ioSlipTime);
  
  // setup StepSlipFn
  StepSlipFn slipfn;
  slipfn.dbFinalSlip(&dbFinalSlip);
  slipfn.dbSlipTime(&dbSlipTime);
  
  spatialdata::units::Nondimensional normalizer;
  const PylithScalar originTime = 5.353;
  
  slipfn.initialize(faultMesh, normalizer, originTime);

  CPPUNIT_ASSERT(slipfn._parameters);
  topology::VecVisitorMesh finalSlipVisitor(slipfn._parameters->get("final slip"));
  const PetscScalar* finalSlipArray = finalSlipVisitor.localArray();CPPUNIT_ASSERT(finalSlipArray);

  topology::VecVisitorMesh slipTimeVisitor(slipfn._parameters->get("slip time"));
  const PetscScalar* slipTimeArray = slipTimeVisitor.localArray();CPPUNIT_ASSERT(slipTimeArray);

  PetscDM faultDMMesh = faultMesh.dmMesh();CPPUNIT_ASSERT(faultDMMesh);
  topology::Stratum verticesStratum(faultDMMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt v = vStart, iPoint = 0; v < vEnd; ++v, ++iPoint) {
    const PetscInt fsoff = finalSlipVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, finalSlipVisitor.sectionDof(v));

    const PetscInt stoff = slipTimeVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(1, slipTimeVisitor.sectionDof(v));

    for(PetscInt d = 0; d < spaceDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.finalSlipE[iPoint*spaceDim+d], finalSlipArray[fsoff+d], tolerance);
    } // for
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime, slipTimeArray[stoff], tolerance);
  } // for

  PYLITH_METHOD_END;
} // _testInitialize


// End of file 
