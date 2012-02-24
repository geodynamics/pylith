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

#include "TestBruneSlipFn.hh" // Implementation of class methods

#include "pylith/faults/BruneSlipFn.hh" // USES BruneSlipFn

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
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestBruneSlipFn );

// ----------------------------------------------------------------------
namespace pylith {
  namespace faults {
    namespace _TestBruneSlipFn {
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
    } // _TestBruneSlipFn
  } // faults
} // pylith

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestBruneSlipFn::testConstructor(void)
{ // testConstructor
  BruneSlipFn slipfn;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbFinalSlip().
void
pylith::faults::TestBruneSlipFn::testDbFinalSlip(void)
{ // testDbFinalSlip
  const char* label = "database ABC";
  BruneSlipFn slipfn;
  
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
pylith::faults::TestBruneSlipFn::testDbSlipTime(void)
{ // testDbSlipTime
  const char* label = "database ABCD";
  BruneSlipFn slipfn;
  
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
pylith::faults::TestBruneSlipFn::testDbRiseTime(void)
{ // testDbRiseTime
  const char* label = "database ABCDE";
  BruneSlipFn slipfn;
  
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
pylith::faults::TestBruneSlipFn::testInitialize1D(void)
{ // testInitialize1D
  const char* meshFilename = "data/line2.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/line2_finalslip.spatialdb";
  const char* slipTimeFilename = "data/line2_sliptime.spatialdb";
  const char* riseTimeFilename = "data/line2_risetime.spatialdb";
  const int constraintPts[] = { 3 };
  const PylithScalar finalSlipE[] = { 2.3 };
  const PylithScalar slipTimeE[] = { 1.2 };
  const PylithScalar riseTimeE[] = { 1.4 };
  const int numConstraintPts = 1;

  _TestBruneSlipFn::DataStruct data = {meshFilename,
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
pylith::faults::TestBruneSlipFn::testInitialize2D(void)
{ // testInitialize2D
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

  _TestBruneSlipFn::DataStruct data = {meshFilename,
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
pylith::faults::TestBruneSlipFn::testInitialize3D(void)
{ // testInitialize3D
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

  _TestBruneSlipFn::DataStruct data = {meshFilename,
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
pylith::faults::TestBruneSlipFn::testSlip(void)
{ // testSlip
  const PylithScalar finalSlipE[] = { 2.3, 0.1, 
				0.0, 0.0};
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const PylithScalar riseTimeE[] = { 1.4, 1.5 };
  const PylithScalar originTime = 5.064;

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  BruneSlipFn slipfn;
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

  const PylithScalar t = 2.134;
  slipfn.slip(&slip, originTime+t);

  const PylithScalar tolerance = 1.0e-06;
  int iPoint = 0;
  const ALE::Obj<RealSection>& slipSection = slip.section();
  CPPUNIT_ASSERT(!slipSection.isNull());
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    PylithScalar slipMag = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      slipMag += pow(finalSlipE[iPoint*spaceDim+iDim], 2);
    slipMag = sqrt(slipMag);
    const PylithScalar peakRate = slipMag / riseTimeE[iPoint] * 1.745;
    const PylithScalar tau = 
      (slipMag > 0.0) ? slipMag / (exp(1.0) * peakRate) : 1.0;
    const PylithScalar t0 = slipTimeE[iPoint];
    const PylithScalar slipNorm = 1.0 - exp(-(t-t0)/tau) * (1.0 + (t-t0)/tau);
    const int fiberDim = slipSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const PylithScalar* vals = slipSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vals);

    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const PylithScalar slipE = finalSlipE[iPoint*spaceDim+iDim] * slipNorm;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, vals[iDim], tolerance);
    } // for
  } // for
} // testSlip

// ----------------------------------------------------------------------
// Test slipIncr().
void
pylith::faults::TestBruneSlipFn::testSlipIncr(void)
{ // testSlipIncr
  const PylithScalar finalSlipE[] = { 2.3, 0.1, 
				0.0, 0.0};
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const PylithScalar riseTimeE[] = { 1.4, 1.5 };
  const PylithScalar originTime = 1.064;

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  BruneSlipFn slipfn;
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

  const PylithScalar t0 = 1.234;
  const PylithScalar t1 = 3.635;
  slipfn.slipIncr(&slip, originTime+t0, originTime+t1);

  const PylithScalar tolerance = 1.0e-06;
  int iPoint = 0;
  const ALE::Obj<RealSection>& slipSection = slip.section();
  CPPUNIT_ASSERT(!slipSection.isNull());
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    PylithScalar slipMag = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      slipMag += pow(finalSlipE[iPoint*spaceDim+iDim], 2);
    slipMag = sqrt(slipMag);
    const PylithScalar peakRate = slipMag / riseTimeE[iPoint] * 1.745;
    const PylithScalar tau = 
      (slipMag > 0.0) ? slipMag / (exp(1.0) * peakRate) : 1.0;
    const PylithScalar tRef = slipTimeE[iPoint];
    const PylithScalar slipNorm0 = 
      (t0 > tRef) ? 1.0 - exp(-(t0-tRef)/tau) * (1.0 + (t0-tRef)/tau) : 0.0;
    const PylithScalar slipNorm1 =
      (t1 > tRef) ? 1.0 - exp(-(t1-tRef)/tau) * (1.0 + (t1-tRef)/tau) : 0.0;

    const int fiberDim = slipSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const PylithScalar* vals = slipSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vals);

    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const PylithScalar slipE = 
	finalSlipE[iPoint*spaceDim+iDim] * (slipNorm1-slipNorm0);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, vals[iDim], tolerance);
    } // for
  } // for
} // testSlipIncr

// ----------------------------------------------------------------------
// Test _slip().
void
pylith::faults::TestBruneSlipFn::testSlipTH(void)
{ // testSlipTH
  const PylithScalar t = 0.734;
  const PylithScalar finalSlip = 4.64;
  const PylithScalar riseTime = 3.23;

  const PylithScalar peakRate = finalSlip / riseTime * 1.745;
  const PylithScalar tau = finalSlip / (exp(1.0) * peakRate);
  const PylithScalar slipE = finalSlip * (1.0 - exp(-t/tau) * (1.0 + t/tau));

  PylithScalar slip = BruneSlipFn::_slipFn(t, finalSlip, riseTime);

  const PylithScalar tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, slip, tolerance);

  slip = BruneSlipFn::_slipFn(-0.5, finalSlip, riseTime);
  CPPUNIT_ASSERT_EQUAL(PylithScalar(0.0), slip);

  slip = BruneSlipFn::_slipFn(1.0e+10, finalSlip, riseTime);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(finalSlip, slip, tolerance);
} // testSlipTH

// ----------------------------------------------------------------------
// Initialize BruneSlipFn.
void
pylith::faults::TestBruneSlipFn::_initialize(topology::Mesh* mesh,
					     topology::SubMesh* faultMesh,
					     BruneSlipFn* slipfn,
					     const PylithScalar originTime)
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != faultMesh);
  CPPUNIT_ASSERT(0 != slipfn);

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
                           faultId,
                           firstFaultVertex, firstLagrangeVertex, firstFaultCell,
                           useLagrangeConstraints);
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

  // setup BruneSlipFn
  slipfn->dbFinalSlip(&dbFinalSlip);
  slipfn->dbSlipTime(&dbSlipTime);
  slipfn->dbRiseTime(&dbRiseTime);
  
  slipfn->initialize(*faultMesh, normalizer, originTime);
} // _initialize

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestBruneSlipFn::_testInitialize(const _TestBruneSlipFn::DataStruct& data)
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
                           data.faultId,
                           firstFaultVertex, firstLagrangeVertex, firstFaultCell,
                           useLagrangeConstraints);
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

  // setup BruneSlipFn
  BruneSlipFn slipfn;
  slipfn.dbFinalSlip(&dbFinalSlip);
  slipfn.dbSlipTime(&dbSlipTime);
  slipfn.dbRiseTime(&dbRiseTime);
  
  spatialdata::units::Nondimensional normalizer;
  const PylithScalar originTime = 5.353;
  
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

  const PylithScalar tolerance = 1.0e-06;
  int iPoint = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    CPPUNIT_ASSERT_EQUAL(spaceDim, finalSlipSection->getFiberDimension(*v_iter));
    const PylithScalar* finalSlipVertex = finalSlipSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != finalSlipVertex);
    for (int iDim=0; iDim < spaceDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.finalSlipE[iPoint*spaceDim+iDim],
				   finalSlipVertex[iDim],
				   tolerance);

    CPPUNIT_ASSERT_EQUAL(1, slipTimeSection->getFiberDimension(*v_iter));
    const PylithScalar* slipTimeVertex = slipTimeSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != slipTimeVertex);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime,
				 slipTimeVertex[0], tolerance);

    CPPUNIT_ASSERT_EQUAL(1, riseTimeSection->getFiberDimension(*v_iter));
    const PylithScalar* riseTimeVertex = riseTimeSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != riseTimeVertex);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.riseTimeE[iPoint],
				 riseTimeVertex[0], tolerance);
  } // for
} // _testInitialize



// End of file 
