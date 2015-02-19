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

#include "TestTimeHistorySlipFn.hh" // Implementation of class methods

#include "pylith/faults/TimeHistorySlipFn.hh" // USES TimeHistorySlipFn

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
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestTimeHistorySlipFn );

// ----------------------------------------------------------------------
namespace pylith {
  namespace faults {
    namespace _TestTimeHistorySlipFn {
      struct DataStruct {
	const char* meshFilename;
	const char* faultLabel;
	const int faultId;
	const char* finalSlipFilename;
	const char* slipTimeFilename;
	const char* timeHistoryFilename;
	const int* constraintPts;
	const PylithScalar* amplitudeE;
	const PylithScalar* slipTimeE;
	const int numConstraintPts;
      }; // DataStruct
    } // _TestTimeHistorySlipFn
  } // faults
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestTimeHistorySlipFn::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  TimeHistorySlipFn slipfn;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbAmplitude().
void
pylith::faults::TestTimeHistorySlipFn::testDbAmplitude(void)
{ // testDbAmplitude
  PYLITH_METHOD_BEGIN;

  const char* label = "database ABC";
  TimeHistorySlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbAmplitude(&db);

  CPPUNIT_ASSERT(slipfn._dbAmplitude);
  CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(slipfn._dbAmplitude->label()));
  CPPUNIT_ASSERT(!slipfn._dbSlipTime);
  CPPUNIT_ASSERT(!slipfn._dbTimeHistory);

  PYLITH_METHOD_END;
} // testDbAmplitude

// ----------------------------------------------------------------------
// Test dbSlipTime().
void
pylith::faults::TestTimeHistorySlipFn::testDbSlipTime(void)
{ // testDbSlipTime
  PYLITH_METHOD_BEGIN;

  const char* label = "database ABCD";
  TimeHistorySlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbSlipTime(&db);

  CPPUNIT_ASSERT(slipfn._dbSlipTime);
  CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(slipfn._dbSlipTime->label()));
  CPPUNIT_ASSERT(!slipfn._dbAmplitude);
  CPPUNIT_ASSERT(!slipfn._dbTimeHistory);

  PYLITH_METHOD_END;
} // testDbSlipTime

// ----------------------------------------------------------------------
// Test dbTimeHistory().
void
pylith::faults::TestTimeHistorySlipFn::testDbTimeHistory(void)
{ // testDbTimeHistory
  PYLITH_METHOD_BEGIN;

  const char* label = "database ABCDE";
  TimeHistorySlipFn slipfn;
  
  spatialdata::spatialdb::TimeHistory db(label);
  slipfn.dbTimeHistory(&db);

  CPPUNIT_ASSERT(slipfn._dbTimeHistory);
  CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(slipfn._dbTimeHistory->label()));
  CPPUNIT_ASSERT(!slipfn._dbAmplitude);
  CPPUNIT_ASSERT(!slipfn._dbSlipTime);

  PYLITH_METHOD_END;
} // testDbTimeHistory

// ----------------------------------------------------------------------
// Test initialize() in 2-D.
void
pylith::faults::TestTimeHistorySlipFn::testInitialize2D(void)
{ // testInitialize2D
  PYLITH_METHOD_BEGIN;

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* timeHistoryFilename = "data/slipfn.timedb";
  const int constraintPts[] = { 3, 4 };
  const PylithScalar amplitudeE[] = { 2.3, 0.1, 
				2.4, 0.2};
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const int numConstraintPts = 2;

  _TestTimeHistorySlipFn::DataStruct data = {meshFilename,
				       faultLabel,
				       faultId,
				       finalSlipFilename,
				       slipTimeFilename,
				       timeHistoryFilename,
				       constraintPts,
				       amplitudeE,
				       slipTimeE,
				       numConstraintPts};
  _testInitialize(data);

  PYLITH_METHOD_END;
} // testInitialize2D

// ----------------------------------------------------------------------
// Test initialize() in 3-D.
void
pylith::faults::TestTimeHistorySlipFn::testInitialize3D(void)
{ // testInitialize3D
  PYLITH_METHOD_BEGIN;

  const char* meshFilename = "data/tet4.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tet4_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tet4_sliptime.spatialdb";
  const char* timeHistoryFilename = "data/slipfn.timedb";
  const int constraintPts[] = { 3, 4, 5 };
  const PylithScalar amplitudeE[] = { 2.3, -0.7, 0.1,
				2.4, -0.8, 0.2,
				2.5, -0.9, 0.3 };
  const PylithScalar slipTimeE[] = { 1.2, 1.3, 1.4 };
  const int numConstraintPts = 3;

  _TestTimeHistorySlipFn::DataStruct data = {meshFilename,
				       faultLabel,
				       faultId,
				       finalSlipFilename,
				       slipTimeFilename,
				       timeHistoryFilename,
				       constraintPts,
				       amplitudeE,
				       slipTimeE,
				       numConstraintPts};
  _testInitialize(data);

  PYLITH_METHOD_END;
} // testInitialize3D

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestTimeHistorySlipFn::testSlip(void)
{ // testSlip
  PYLITH_METHOD_BEGIN;

  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const PylithScalar slipE[] = { 0.92, 0.04,
			   0.84, 0.07 };
  const PylithScalar originTime = 5.064;

  topology::Mesh mesh;
  topology::Mesh faultMesh;
  TimeHistorySlipFn slipfn;
  spatialdata::spatialdb::TimeHistory th;
  _initialize(&mesh, &faultMesh, &slipfn, &th, originTime);
  
  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();CPPUNIT_ASSERT(cs);
  const int spaceDim = cs->spaceDim();

  topology::Field slip(faultMesh);
  slip.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  slip.allocate();

  const PylithScalar t = 2.0;
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
    
    for(PetscInt d = 0; d < spaceDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE[iPoint*spaceDim+d], slipArray[off+d], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testSlip

// ----------------------------------------------------------------------
// Initialize TimeHistorySlipFn.
void
pylith::faults::TestTimeHistorySlipFn::_initialize(topology::Mesh* mesh,
						   topology::Mesh* faultMesh,
						   TimeHistorySlipFn* slipfn,
						   spatialdata::spatialdb::TimeHistory* th,
						   const PylithScalar originTime)
{ // _initialize
  PYLITH_METHOD_BEGIN;

  assert(slipfn);
  PetscErrorCode err;

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* timeHistoryFilename = "data/slipfn.timedb";

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
  spatialdata::spatialdb::SimpleDB dbAmplitude("slip amplitude");
  spatialdata::spatialdb::SimpleIOAscii ioFinalSlip;
  ioFinalSlip.filename(finalSlipFilename);
  dbAmplitude.ioHandler(&ioFinalSlip);
  
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
  ioSlipTime.filename(slipTimeFilename);
  dbSlipTime.ioHandler(&ioSlipTime);
  
  th->label("time history");
  th->filename(timeHistoryFilename);

  spatialdata::units::Nondimensional normalizer;

  // setup TimeHistorySlipFn
  slipfn->dbAmplitude(&dbAmplitude);
  slipfn->dbSlipTime(&dbSlipTime);
  slipfn->dbTimeHistory(th);
  
  slipfn->initialize(*faultMesh, normalizer, originTime);

  PYLITH_METHOD_END;
} // _initialize

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestTimeHistorySlipFn::_testInitialize(const _TestTimeHistorySlipFn::DataStruct& data)
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
  spatialdata::spatialdb::SimpleDB dbAmplitude("slip amplitude");
  spatialdata::spatialdb::SimpleIOAscii ioFinalSlip;
  ioFinalSlip.filename(data.finalSlipFilename);
  dbAmplitude.ioHandler(&ioFinalSlip);
  
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
  ioSlipTime.filename(data.slipTimeFilename);
  dbSlipTime.ioHandler(&ioSlipTime);
  
  spatialdata::spatialdb::TimeHistory dbTimeHistory("time history");
  dbTimeHistory.filename(data.timeHistoryFilename);

  // setup TimeHistorySlipFn
  TimeHistorySlipFn slipfn;
  slipfn.dbAmplitude(&dbAmplitude);
  slipfn.dbSlipTime(&dbSlipTime);
  slipfn.dbTimeHistory(&dbTimeHistory);
  
  spatialdata::units::Nondimensional normalizer;
  const PylithScalar originTime = 5.353;
  
  slipfn.initialize(faultMesh, normalizer, originTime);

  CPPUNIT_ASSERT(slipfn._parameters);
  topology::VecVisitorMesh finalSlipVisitor(slipfn._parameters->get("slip amplitude"));
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
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.amplitudeE[iPoint*spaceDim+d], finalSlipArray[fsoff+d], tolerance);
    } // for
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime, slipTimeArray[stoff], tolerance);
  } // for

  PYLITH_METHOD_END;
} // _testInitialize


// End of file 
