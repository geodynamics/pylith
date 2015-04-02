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

#include "TestLiuCosSlipFn.hh" // Implementation of class methods

#include "pylith/faults/LiuCosSlipFn.hh" // USES LiuCosSlipFn

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
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestLiuCosSlipFn );

// ----------------------------------------------------------------------
namespace pylith {
  namespace faults {
    namespace _TestLiuCosSlipFn {
      struct DataStruct {
	const char* meshFilename;
	const char* faultLabel;
	const int faultId;
	const char* finalSlipFilename;
	const char* slipTimeFilename;
	const char* riseTimeFilename;
	const int* constraintPts;
	const PylithScalar* finalSlipE;
	const PylithScalar* slipTimeE;
	const PylithScalar* riseTimeE;
	const int numConstraintPts;
      }; // DataStruct
    } // _TestLiuCosSlipFn
  } // faults
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestLiuCosSlipFn::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  LiuCosSlipFn slipfn;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbFinalSlip().
void
pylith::faults::TestLiuCosSlipFn::testDbFinalSlip(void)
{ // testDbFinalSlip
  PYLITH_METHOD_BEGIN;

  const char* label = "database ABC";
  LiuCosSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbFinalSlip(&db);

  CPPUNIT_ASSERT(slipfn._dbFinalSlip);
  CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(slipfn._dbFinalSlip->label()));
  CPPUNIT_ASSERT(!slipfn._dbSlipTime);
  CPPUNIT_ASSERT(!slipfn._dbRiseTime);

  PYLITH_METHOD_END;
} // testDbFinalSlip

// ----------------------------------------------------------------------
// Test dbSlipTime().
void
pylith::faults::TestLiuCosSlipFn::testDbSlipTime(void)
{ // testDbSlipTime
  PYLITH_METHOD_BEGIN;

  const char* label = "database ABCD";
  LiuCosSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbSlipTime(&db);

  CPPUNIT_ASSERT(slipfn._dbSlipTime);
  CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(slipfn._dbSlipTime->label()));
  CPPUNIT_ASSERT(!slipfn._dbFinalSlip);
  CPPUNIT_ASSERT(!slipfn._dbRiseTime);

  PYLITH_METHOD_END;
} // testDbSlipTime

// ----------------------------------------------------------------------
// Test dbRiseTime().
void
pylith::faults::TestLiuCosSlipFn::testDbRiseTime(void)
{ // testDbRiseTime
  PYLITH_METHOD_BEGIN;

  const char* label = "database ABCDE";
  LiuCosSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbRiseTime(&db);

  CPPUNIT_ASSERT(slipfn._dbRiseTime);
  CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(slipfn._dbRiseTime->label()));
  CPPUNIT_ASSERT(!slipfn._dbFinalSlip);
  CPPUNIT_ASSERT(!slipfn._dbSlipTime);

  PYLITH_METHOD_END;
} // testDbRiseTime

// ----------------------------------------------------------------------
// Test initialize() in 2-D.
void
pylith::faults::TestLiuCosSlipFn::testInitialize2D(void)
{ // testInitialize2D
  PYLITH_METHOD_BEGIN;

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* riseTimeFilename = "data/tri3_risetime.spatialdb";
  const int constraintPts[] = { 3, 4 };
  const PylithScalar finalSlipE[] = { 2.3, 0.1, 
				2.4, 0.2};
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const PylithScalar riseTimeE[] = { 1.4, 1.5 };
  const int numConstraintPts = 2;

  _TestLiuCosSlipFn::DataStruct data = {meshFilename,
				       faultLabel,
				       faultId,
				       finalSlipFilename,
				       slipTimeFilename,
				       riseTimeFilename,
				       constraintPts,
				       finalSlipE,
				       slipTimeE,
				       riseTimeE,
				       numConstraintPts};
  _testInitialize(data);

  PYLITH_METHOD_END;
} // testInitialize2D

// ----------------------------------------------------------------------
// Test initialize() in 3-D.
void
pylith::faults::TestLiuCosSlipFn::testInitialize3D(void)
{ // testInitialize3D
  PYLITH_METHOD_BEGIN;

  const char* meshFilename = "data/tet4.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tet4_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tet4_sliptime.spatialdb";
  const char* riseTimeFilename = "data/tet4_risetime.spatialdb";
  const int constraintPts[] = { 3, 4, 5 };
  const PylithScalar finalSlipE[] = { 2.3, -0.7, 0.1,
				2.4, -0.8, 0.2,
				2.5, -0.9, 0.3 };
  const PylithScalar slipTimeE[] = { 1.2, 1.3, 1.4 };
  const PylithScalar riseTimeE[] = { 1.5, 1.6, 1.7 };
  const int numConstraintPts = 3;

  _TestLiuCosSlipFn::DataStruct data = {meshFilename,
				       faultLabel,
				       faultId,
				       finalSlipFilename,
				       slipTimeFilename,
				       riseTimeFilename,
				       constraintPts,
				       finalSlipE,
				       slipTimeE,
				       riseTimeE,
				       numConstraintPts};
  _testInitialize(data);

  PYLITH_METHOD_END;
} // testInitialize3D

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestLiuCosSlipFn::testSlip(void)
{ // testSlip
  const PylithScalar finalSlipE[] = { 2.3, 0.1, 
				0.0, 0.0};
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const PylithScalar riseTimeE[] = { 1.4, 1.5 };
  const PylithScalar originTime = 5.064;

  topology::Mesh mesh;
  topology::Mesh faultMesh;
  LiuCosSlipFn slipfn;
  _initialize(&mesh, &faultMesh, &slipfn, originTime);
  
  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();CPPUNIT_ASSERT(cs);
  const int spaceDim = cs->spaceDim();

  topology::Field slip(faultMesh);
  slip.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  slip.allocate();

  const PylithScalar t = 2.134;
  slipfn.slip(&slip, originTime+t);

  PetscDM dmMesh = faultMesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  topology::VecVisitorMesh slipVisitor(slip);
  const PetscScalar* slipArray = slipVisitor.localArray();CPPUNIT_ASSERT(slipArray);

  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt v = vStart, iPoint = 0; v < vEnd; ++v, ++iPoint) {
    PylithScalar slipMag = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      slipMag += pow(finalSlipE[iPoint*spaceDim+iDim], 2);
    } // for
    slipMag = sqrt(slipMag);

    const PylithScalar slipNorm = (slipMag > 0.0) ?
      _slipFn(t - slipTimeE[iPoint], slipMag, riseTimeE[iPoint]) / slipMag : 0.0;

    const PetscInt off = slipVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, slipVisitor.sectionDof(v));

    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar slipE = finalSlipE[iPoint*spaceDim+d] * slipNorm;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, slipArray[off+d], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testSlip

// ----------------------------------------------------------------------
// Test _slip().
void
pylith::faults::TestLiuCosSlipFn::testSlipTH(void)
{ // testSlipTH
  PYLITH_METHOD_BEGIN;

  const PylithScalar t = 0.734;
  const PylithScalar finalSlip = 4.64;
  const PylithScalar riseTime = 3.23;

  const PylithScalar slipE = _slipFn(t, finalSlip, riseTime);

  PylithScalar slip = LiuCosSlipFn::_slipFn(t, finalSlip, riseTime);

  const PylithScalar tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, slip, tolerance);

  slip = LiuCosSlipFn::_slipFn(-0.5, finalSlip, riseTime);
  CPPUNIT_ASSERT_EQUAL(PylithScalar(0.0), slip);

  slip = LiuCosSlipFn::_slipFn(1.0e+10, finalSlip, riseTime);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(finalSlip, slip, tolerance);

  PYLITH_METHOD_END;
} // testSlipTH

// ----------------------------------------------------------------------
// Initialize LiuCosSlipFn.
void
pylith::faults::TestLiuCosSlipFn::_initialize(topology::Mesh* mesh,
					      topology::Mesh* faultMesh,
					      LiuCosSlipFn* slipfn,
					      const PylithScalar originTime)
{ // _initialize
  PYLITH_METHOD_BEGIN;

  assert(slipfn);
  PetscErrorCode  err;

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslipB.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* riseTimeFilename = "data/tri3_risetime.spatialdb";

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
  
  spatialdata::spatialdb::SimpleDB dbRiseTime("rise time");
  spatialdata::spatialdb::SimpleIOAscii ioRiseTime;
  ioRiseTime.filename(riseTimeFilename);
  dbRiseTime.ioHandler(&ioRiseTime);

  spatialdata::units::Nondimensional normalizer;

  // setup LiuCosSlipFn
  slipfn->dbFinalSlip(&dbFinalSlip);
  slipfn->dbSlipTime(&dbSlipTime);
  slipfn->dbRiseTime(&dbRiseTime);
  
  slipfn->initialize(*faultMesh, normalizer, originTime);

  PYLITH_METHOD_END;
} // _initialize

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestLiuCosSlipFn::_testInitialize(const _TestLiuCosSlipFn::DataStruct& data)
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
  
  spatialdata::spatialdb::SimpleDB dbRiseTime("rise time");
  spatialdata::spatialdb::SimpleIOAscii ioRiseTime;
  ioRiseTime.filename(data.riseTimeFilename);
  dbRiseTime.ioHandler(&ioRiseTime);

  // setup LiuCosSlipFn
  LiuCosSlipFn slipfn;
  slipfn.dbFinalSlip(&dbFinalSlip);
  slipfn.dbSlipTime(&dbSlipTime);
  slipfn.dbRiseTime(&dbRiseTime);
  
  spatialdata::units::Nondimensional normalizer;
  const PylithScalar originTime = 5.353;
  
  slipfn.initialize(faultMesh, normalizer, originTime);

  CPPUNIT_ASSERT(slipfn._parameters);
  topology::VecVisitorMesh finalSlipVisitor(slipfn._parameters->get("final slip"));
  const PetscScalar* finalSlipArray = finalSlipVisitor.localArray();CPPUNIT_ASSERT(finalSlipArray);

  topology::VecVisitorMesh slipTimeVisitor(slipfn._parameters->get("slip time"));
  const PetscScalar* slipTimeArray = slipTimeVisitor.localArray();CPPUNIT_ASSERT(slipTimeArray);

  topology::VecVisitorMesh riseTimeVisitor(slipfn._parameters->get("rise time"));
  const PetscScalar* riseTimeArray = riseTimeVisitor.localArray();CPPUNIT_ASSERT(riseTimeArray);

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

    const PetscInt rtoff = riseTimeVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(1, riseTimeVisitor.sectionDof(v));

    for(PetscInt d = 0; d < spaceDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.finalSlipE[iPoint*spaceDim+d], finalSlipArray[fsoff+d], tolerance);
    } // for
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime, slipTimeArray[stoff], tolerance);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.riseTimeE[iPoint], riseTimeArray[rtoff], tolerance);
  } // for

  PYLITH_METHOD_END;
} // _testInitialize

// ----------------------------------------------------------------------
// Slip time function.
PylithScalar
pylith::faults::TestLiuCosSlipFn::_slipFn(const PylithScalar t,
					  const PylithScalar finalSlip,
					  const PylithScalar riseTime)
{ // _slipFn
  PYLITH_METHOD_BEGIN;

  const float tau = riseTime * 1.525;
  const float tau1 = 0.13 * tau;
  const float tau2 = tau - tau1;
  const float Cn = 
    M_PI /  (1.4 * M_PI * tau1 + 1.2 * tau1 + 0.3 * M_PI * tau2);
  
  PylithScalar slip = 0.0;
  if (t <= tau1) {
    slip = 0.7*t - 0.7*tau1/M_PI*sin(M_PI*t/tau1) 
      - 0.6*tau1/(0.5*M_PI)*(cos(0.5*M_PI*t/tau1) - 1.0);
    slip *= Cn;
  } else if (t <= 2.0*tau1) {
    slip = 1.0*t - 0.7*tau1/M_PI*sin(M_PI*t/tau1) + 
      0.3*tau2/M_PI*sin(M_PI*(t-tau1)/tau2) + 
      1.2*tau1/M_PI - 0.3*tau1;
    slip *= Cn;
  } else if (t <= tau) {
    slip = 0.3*t + 0.3*tau2/M_PI*sin(M_PI*(t-tau1)/tau2) + 
      1.1*tau1 + 1.2*tau1/M_PI;
    slip *= Cn;	
  } else
    slip = 1.0;
  slip *= finalSlip;

  PYLITH_METHOD_RETURN(slip);
} // _slipFn


// End of file 
