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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestLiuCosSlipFn.hh" // Implementation of class methods

#include "pylith/faults/LiuCosSlipFn.hh" // USES LiuCosSlipFn

#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestLiuCosSlipFn );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

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
	const double* finalSlipE;
	const double* slipTimeE;
	const double* riseTimeE;
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
  LiuCosSlipFn slipfn;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbFinalSlip().
void
pylith::faults::TestLiuCosSlipFn::testDbFinalSlip(void)
{ // testDbFinalSlip
  const char* label = "database ABC";
  LiuCosSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbFinalSlip(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbFinalSlip);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbFinalSlip->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbSlipTime);
  CPPUNIT_ASSERT(0 == slipfn._dbRiseTime);
} // testDbFinalSlip

// ----------------------------------------------------------------------
// Test dbSlipTime().
void
pylith::faults::TestLiuCosSlipFn::testDbSlipTime(void)
{ // testDbSlipTime
  const char* label = "database ABCD";
  LiuCosSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbSlipTime(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbSlipTime);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbSlipTime->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbFinalSlip);
  CPPUNIT_ASSERT(0 == slipfn._dbRiseTime);
} // testDbSlipTime

// ----------------------------------------------------------------------
// Test dbRiseTime().
void
pylith::faults::TestLiuCosSlipFn::testDbRiseTime(void)
{ // testDbRiseTime
  const char* label = "database ABCDE";
  LiuCosSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbRiseTime(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbRiseTime);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbRiseTime->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbFinalSlip);
  CPPUNIT_ASSERT(0 == slipfn._dbSlipTime);
} // testDbRiseTime

// ----------------------------------------------------------------------
// Test initialize() in 1-D.
void
pylith::faults::TestLiuCosSlipFn::testInitialize1D(void)
{ // testInitialize1D
  const char* meshFilename = "data/line2.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/line2_finalslip.spatialdb";
  const char* slipTimeFilename = "data/line2_sliptime.spatialdb";
  const char* riseTimeFilename = "data/line2_risetime.spatialdb";
  const int constraintPts[] = { 3 };
  const double finalSlipE[] = { 2.3 };
  const double slipTimeE[] = { 1.2 };
  const double riseTimeE[] = { 1.4 };
  const int numConstraintPts = 1;

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
} // testInitialize1D

// ----------------------------------------------------------------------
// Test initialize() in 2-D.
void
pylith::faults::TestLiuCosSlipFn::testInitialize2D(void)
{ // testInitialize2D
  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* riseTimeFilename = "data/tri3_risetime.spatialdb";
  const int constraintPts[] = { 3, 4 };
  const double finalSlipE[] = { 2.3, 0.1, 
				2.4, 0.2};
  const double slipTimeE[] = { 1.2, 1.3 };
  const double riseTimeE[] = { 1.4, 1.5 };
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
} // testInitialize2D

// ----------------------------------------------------------------------
// Test initialize() in 3-D.
void
pylith::faults::TestLiuCosSlipFn::testInitialize3D(void)
{ // testInitialize3D
  const char* meshFilename = "data/tet4.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tet4_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tet4_sliptime.spatialdb";
  const char* riseTimeFilename = "data/tet4_risetime.spatialdb";
  const int constraintPts[] = { 3, 4, 5 };
  const double finalSlipE[] = { 2.3, -0.7, 0.1,
				2.4, -0.8, 0.2,
				2.5, -0.9, 0.3 };
  const double slipTimeE[] = { 1.2, 1.3, 1.4 };
  const double riseTimeE[] = { 1.5, 1.6, 1.7 };
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
} // testInitialize3D

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestLiuCosSlipFn::testSlip(void)
{ // testSlip
  const double finalSlipE[] = { 2.3, 0.1, 
				0.0, 0.0};
  const double slipTimeE[] = { 1.2, 1.3 };
  const double riseTimeE[] = { 1.4, 1.5 };
  const double originTime = 5.064;

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  LiuCosSlipFn slipfn;
  _initialize(&mesh, &faultMesh, &slipfn, originTime);
  
  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  CPPUNIT_ASSERT(0 != cs);

  const int spaceDim = cs->spaceDim();
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh.sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
    faultSieveMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  topology::Field<topology::SubMesh> slip(faultMesh);
  slip.newSection(vertices, spaceDim);
  slip.allocate();

  const double t = 2.134;
  slipfn.slip(&slip, originTime+t);

  const double tolerance = 1.0e-06;
  int iPoint = 0;
  const ALE::Obj<RealSection>& slipSection = slip.section();
  CPPUNIT_ASSERT(!slipSection.isNull());
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    double slipMag = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      slipMag += pow(finalSlipE[iPoint*spaceDim+iDim], 2);
    slipMag = sqrt(slipMag);

    const double slipNorm = (slipMag > 0.0) ?
      _slipFn(t - slipTimeE[iPoint], slipMag, riseTimeE[iPoint]) / slipMag :
      0.0;

    const int fiberDim = slipSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const double* vals = slipSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vals);
    
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const double slipE = finalSlipE[iPoint*spaceDim+iDim] * slipNorm;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, vals[iDim], tolerance);
    } // for
  } // for
} // testSlip

// ----------------------------------------------------------------------
// Test slipIncr().
void
pylith::faults::TestLiuCosSlipFn::testSlipIncr(void)
{ // testSlipIncr
  const double finalSlipE[] = { 2.3, 0.1, 
				0.0, 0.0};
  const double slipTimeE[] = { 1.2, 1.3 };
  const double riseTimeE[] = { 1.4, 1.5 };
  const double originTime = 1.064;

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  LiuCosSlipFn slipfn;
  _initialize(&mesh, &faultMesh, &slipfn, originTime);

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  CPPUNIT_ASSERT(0 != cs);

  const int spaceDim = cs->spaceDim();
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh.sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
    faultSieveMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  topology::Field<topology::SubMesh> slip(faultMesh);
  slip.newSection(vertices, spaceDim);
  slip.allocate();

  const double t0 = 1.234;
  const double t1 = 3.635;
  slipfn.slipIncr(&slip, originTime+t0, originTime+t1);

  const double tolerance = 1.0e-06;
  int iPoint = 0;
  const ALE::Obj<RealSection>& slipSection = slip.section();
  CPPUNIT_ASSERT(!slipSection.isNull());
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    double slipMag = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      slipMag += pow(finalSlipE[iPoint*spaceDim+iDim], 2);
    slipMag = sqrt(slipMag);

    const double slipNorm0 = (slipMag > 0.0) ?
      _slipFn(t0 - slipTimeE[iPoint], slipMag, riseTimeE[iPoint]) / slipMag :
      0.0;
    const double slipNorm1 = (slipMag > 0.0) ?
      _slipFn(t1 - slipTimeE[iPoint], slipMag, riseTimeE[iPoint]) / slipMag :
      0.0;

    const int fiberDim = slipSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const double* vals = slipSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vals);

    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const double slipE = 
	finalSlipE[iPoint*spaceDim+iDim] * (slipNorm1-slipNorm0);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, vals[iDim], tolerance);
    } // for
  } // for
} // testSlipIncr

// ----------------------------------------------------------------------
// Test _slip().
void
pylith::faults::TestLiuCosSlipFn::testSlipTH(void)
{ // testSlipTH
  const double t = 0.734;
  const double finalSlip = 4.64;
  const double riseTime = 3.23;

  const double slipE = _slipFn(t, finalSlip, riseTime);

  double slip = LiuCosSlipFn::_slipFn(t, finalSlip, riseTime);

  const double tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, slip, tolerance);

  slip = LiuCosSlipFn::_slipFn(-0.5, finalSlip, riseTime);
  CPPUNIT_ASSERT_EQUAL(0.0, slip);

  slip = LiuCosSlipFn::_slipFn(1.0e+10, finalSlip, riseTime);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(finalSlip, slip, tolerance);
} // testSlipTH

// ----------------------------------------------------------------------
// Initialize LiuCosSlipFn.
void
pylith::faults::TestLiuCosSlipFn::_initialize(topology::Mesh* mesh,
					      topology::SubMesh* faultMesh,
					      LiuCosSlipFn* slipfn,
					      const double originTime)
{ // _initialize
  assert(0 != slipfn);

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
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);

  // Create fault mesh
  int firstFaultVertex    = 0;
  int firstLagrangeVertex = mesh->sieveMesh()->getIntSection(faultLabel)->size();
  int firstFaultCell      = mesh->sieveMesh()->getIntSection(faultLabel)->size();
  const bool useLagrangeConstraints = true;
  if (useLagrangeConstraints) {
    firstFaultCell += mesh->sieveMesh()->getIntSection(faultLabel)->size();
  }
  ALE::Obj<SieveFlexMesh> faultBoundary = 0;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CohesiveTopology::createFault(faultMesh, faultBoundary,
                                *mesh, sieveMesh->getIntSection(faultLabel));
  CohesiveTopology::create(mesh, *faultMesh, faultBoundary, 
                           sieveMesh->getIntSection(faultLabel),
                           faultId, firstFaultVertex, firstLagrangeVertex,
			   firstFaultCell, useLagrangeConstraints);
  // Need to copy coordinates from mesh to fault mesh since we are not
  // using create() instead of createParallel().
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh->sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  faultSieveMesh->setRealSection("coordinates", 
				 sieveMesh->getRealSection("coordinates"));

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
} // _initialize

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestLiuCosSlipFn::_testInitialize(const _TestLiuCosSlipFn::DataStruct& data)
{ // _testInitialize
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
  topology::SubMesh faultMesh;
  int firstFaultVertex    = 0;
  int firstLagrangeVertex = mesh.sieveMesh()->getIntSection(data.faultLabel)->size();
  int firstFaultCell      = mesh.sieveMesh()->getIntSection(data.faultLabel)->size();
  const bool useLagrangeConstraints = true;
  if (useLagrangeConstraints) {
    firstFaultCell += mesh.sieveMesh()->getIntSection(data.faultLabel)->size();
  }
  ALE::Obj<SieveFlexMesh> faultBoundary = 0;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CohesiveTopology::createFault(&faultMesh, faultBoundary,
                                mesh, sieveMesh->getIntSection(data.faultLabel));
  CohesiveTopology::create(&mesh, faultMesh, faultBoundary, 
                           sieveMesh->getIntSection(data.faultLabel),
                           data.faultId, firstFaultVertex, firstLagrangeVertex,
			   firstFaultCell, useLagrangeConstraints);
  // Need to copy coordinates from mesh to fault mesh since we are not
  // using create() instead of createParallel().
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh.sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  faultSieveMesh->setRealSection("coordinates", 
				 sieveMesh->getRealSection("coordinates"));

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
  const double originTime = 5.353;
  
  slipfn.initialize(faultMesh, normalizer, originTime);

  const ALE::Obj<SieveMesh::label_sequence>& vertices =
    faultSieveMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

  CPPUNIT_ASSERT(0 != slipfn._parameters);
  const ALE::Obj<RealSection>& finalSlipSection =
    slipfn._parameters->get("final slip").section();
  CPPUNIT_ASSERT(!finalSlipSection.isNull());
  const ALE::Obj<RealSection>& slipTimeSection =
    slipfn._parameters->get("slip time").section();
  CPPUNIT_ASSERT(!slipTimeSection.isNull());
  const ALE::Obj<RealSection>& riseTimeSection =
    slipfn._parameters->get("rise time").section();
  CPPUNIT_ASSERT(!riseTimeSection.isNull());

  const double tolerance = 1.0e-06;
  int iPoint = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    CPPUNIT_ASSERT_EQUAL(spaceDim, finalSlipSection->getFiberDimension(*v_iter));
    const double* finalSlipVertex = finalSlipSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != finalSlipVertex);
    for (int iDim=0; iDim < spaceDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.finalSlipE[iPoint*spaceDim+iDim],
				   finalSlipVertex[iDim],
				   tolerance);

    CPPUNIT_ASSERT_EQUAL(1, slipTimeSection->getFiberDimension(*v_iter));
    const double* slipTimeVertex = slipTimeSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != slipTimeVertex);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime,
				 slipTimeVertex[0], tolerance);

    CPPUNIT_ASSERT_EQUAL(1, riseTimeSection->getFiberDimension(*v_iter));
    const double* riseTimeVertex = riseTimeSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != riseTimeVertex);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.riseTimeE[iPoint],
				 riseTimeVertex[0], tolerance);
  } // for
} // _testInitialize

// ----------------------------------------------------------------------
// Slip time function.
double
pylith::faults::TestLiuCosSlipFn::_slipFn(const double t,
					  const double finalSlip,
					  const double riseTime)
{ // _slipFn
  const float tau = riseTime * 1.525;
  const float tau1 = 0.13 * tau;
  const float tau2 = tau - tau1;
  const float Cn = 
    M_PI /  (1.4 * M_PI * tau1 + 1.2 * tau1 + 0.3 * M_PI * tau2);
  
  double slip = 0.0;
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

  return slip;
} // _slipFn


// End of file 
