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

#include "TestConstRateSlipFn.hh" // Implementation of class methods

#include "pylith/faults/ConstRateSlipFn.hh" // USES ConstRateSlipFn

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
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestConstRateSlipFn );

// ----------------------------------------------------------------------
namespace pylith {
  namespace faults {
    namespace _TestConstRateSlipFn {
      struct DataStruct {
	const char* meshFilename;
	const char* faultLabel;
	const int faultId;
	const char* slipRateFilename;
	const char* slipTimeFilename;
	const int* constraintPts;
	const PylithScalar* slipRateE;
	const PylithScalar* slipTimeE;
	const int numConstraintPts;
      }; // DataStruct
    } // _TestConstRateSlipFn
  } // faults
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestConstRateSlipFn::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  ConstRateSlipFn slipfn;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbSlipRate().
void
pylith::faults::TestConstRateSlipFn::testDbSlipRate(void)
{ // testDbSlipRate
  PYLITH_METHOD_BEGIN;

  const char* label = "database ABC";
  ConstRateSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbSlipRate(&db);

  CPPUNIT_ASSERT(slipfn._dbSlipRate);
  CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(slipfn._dbSlipRate->label()));
  CPPUNIT_ASSERT(!slipfn._dbSlipTime);

  PYLITH_METHOD_END;
} // testDbSlipRate

// ----------------------------------------------------------------------
// Test dbSlipTime().
void
pylith::faults::TestConstRateSlipFn::testDbSlipTime(void)
{ // testDbSlipTime
  PYLITH_METHOD_BEGIN;

 const char* label = "database ABCD";
  ConstRateSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbSlipTime(&db);

  CPPUNIT_ASSERT(slipfn._dbSlipTime);
  CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(slipfn._dbSlipTime->label()));
  CPPUNIT_ASSERT(!slipfn._dbSlipRate);

  PYLITH_METHOD_END;
} // testDbSlipTime

// ----------------------------------------------------------------------
// Test initialize() in 2-D.
void
pylith::faults::TestConstRateSlipFn::testInitialize2D(void)
{ // testInitialize2D
  PYLITH_METHOD_BEGIN;

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* slipRateFilename = "data/tri3_sliprate.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const int constraintPts[] = { 3, 4 };
  const PylithScalar slipRateE[] = { 0.1, 0.2, 
			       0.3, 0.4 };
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const int numConstraintPts = 2;

  _TestConstRateSlipFn::DataStruct data = {meshFilename,
					   faultLabel,
					   faultId,
					   slipRateFilename,
					   slipTimeFilename,
					   constraintPts,
					   slipRateE,
					   slipTimeE,
					   numConstraintPts};
  _testInitialize(data);

  PYLITH_METHOD_END;
} // testInitialize2D

// ----------------------------------------------------------------------
// Test initialize() in 3-D.
void
pylith::faults::TestConstRateSlipFn::testInitialize3D(void)
{ // testInitialize3D
  PYLITH_METHOD_BEGIN;

  const char* meshFilename = "data/tet4.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* slipRateFilename = "data/tet4_sliprate.spatialdb";
  const char* slipTimeFilename = "data/tet4_sliptime.spatialdb";
  const int constraintPts[] = { 3, 4, 5 };
  const PylithScalar slipRateE[] = { 1.6, -0.7, 0.1,
			       1.7, -0.8, 0.2,
			       1.8, -0.9, 0.3 };
  const PylithScalar slipTimeE[] = { 1.2, 1.3, 1.4 };
  const int numConstraintPts = 3;

  _TestConstRateSlipFn::DataStruct data = {meshFilename,
					   faultLabel,
					   faultId,
					   slipRateFilename,
					   slipTimeFilename,
					   constraintPts,
					   slipRateE,
					   slipTimeE,
					   numConstraintPts};
  _testInitialize(data);

  PYLITH_METHOD_END;
} // testInitialize3D

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestConstRateSlipFn::testSlip(void)
{ // testSlip
  PYLITH_METHOD_BEGIN;

  const PylithScalar slipRateE[] = { 0.1, 0.2, 
			       0.3, 0.4};
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const PylithScalar originTime = 5.064;

  topology::Mesh mesh;
  topology::Mesh faultMesh;
  ConstRateSlipFn slipfn;
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
    const PylithScalar t0 = slipTimeE[iPoint];

    const PetscInt off = slipVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, slipVisitor.sectionDof(v));

    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar slipE = (t - t0) > 0.0 ? slipRateE[iPoint*spaceDim+d] * (t - t0) : 0.0;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, slipArray[off+d], tolerance);
    }
  } // for

  PYLITH_METHOD_END;
} // testSlip

// ----------------------------------------------------------------------
// Initialize ConstRateSlipFn.
void
pylith::faults::TestConstRateSlipFn::_initialize(topology::Mesh* mesh,
						 topology::Mesh* faultMesh,
					    	 ConstRateSlipFn* slipfn,
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
  const char* slipRateFilename = "data/tri3_sliprate.spatialdb";
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
  spatialdata::spatialdb::SimpleDB dbSlipRate("slip rate");
  spatialdata::spatialdb::SimpleIOAscii ioSlipRate;
  ioSlipRate.filename(slipRateFilename);
  dbSlipRate.ioHandler(&ioSlipRate);
  
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
  ioSlipTime.filename(slipTimeFilename);
  dbSlipTime.ioHandler(&ioSlipTime);
  
  spatialdata::units::Nondimensional normalizer;

  // setup ConstRateSlipFn
  slipfn->dbSlipRate(&dbSlipRate);
  slipfn->dbSlipTime(&dbSlipTime);
  
  slipfn->initialize(*faultMesh, normalizer, originTime);

  PYLITH_METHOD_END;
} // _initialize

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestConstRateSlipFn::_testInitialize(const _TestConstRateSlipFn::DataStruct& data)
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
  spatialdata::spatialdb::SimpleDB dbSlipRate("slip rate");
  spatialdata::spatialdb::SimpleIOAscii ioSlipRate;
  ioSlipRate.filename(data.slipRateFilename);
  dbSlipRate.ioHandler(&ioSlipRate);
  
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
  ioSlipTime.filename(data.slipTimeFilename);
  dbSlipTime.ioHandler(&ioSlipTime);
  
  // setup ConstRateSlipFn
  ConstRateSlipFn slipfn;
  slipfn.dbSlipRate(&dbSlipRate);
  slipfn.dbSlipTime(&dbSlipTime);
  
  spatialdata::units::Nondimensional normalizer;
  const PylithScalar originTime = 5.353;
  
  slipfn.initialize(faultMesh, normalizer, originTime);

  CPPUNIT_ASSERT(slipfn._parameters);
  topology::VecVisitorMesh slipRateVisitor(slipfn._parameters->get("slip rate"));
  const PetscScalar* slipRateArray = slipRateVisitor.localArray();CPPUNIT_ASSERT(slipRateArray);

  topology::VecVisitorMesh slipTimeVisitor(slipfn._parameters->get("slip time"));
  const PetscScalar* slipTimeArray = slipTimeVisitor.localArray();CPPUNIT_ASSERT(slipTimeArray);

  PetscDM faultDMMesh = faultMesh.dmMesh();CPPUNIT_ASSERT(faultDMMesh);
  topology::Stratum verticesStratum(faultDMMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt v = vStart, iPoint = 0; v < vEnd; ++v, ++iPoint) {
    const PetscInt sroff = slipRateVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, slipRateVisitor.sectionDof(v));

    const PetscInt stoff = slipTimeVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(1, slipTimeVisitor.sectionDof(v));

    for(PetscInt d = 0; d < spaceDim; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipRateE[iPoint*spaceDim+d], slipRateArray[sroff+d], tolerance);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime, slipTimeArray[stoff], tolerance);
  } // for

  PYLITH_METHOD_END;
} // _testInitialize


// End of file 
