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

#include "TestTimeHistorySlipFn.hh" // Implementation of class methods

#include "pylith/faults/TimeHistorySlipFn.hh" // USES TimeHistorySlipFn

#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestTimeHistorySlipFn );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

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
	const double* amplitudeE;
	const double* slipTimeE;
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
  TimeHistorySlipFn slipfn;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbAmplitude().
void
pylith::faults::TestTimeHistorySlipFn::testDbAmplitude(void)
{ // testDbAmplitude
  const char* label = "database ABC";
  TimeHistorySlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbAmplitude(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbAmplitude);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbAmplitude->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbSlipTime);
  CPPUNIT_ASSERT(0 == slipfn._dbTimeHistory);
} // testDbAmplitude

// ----------------------------------------------------------------------
// Test dbSlipTime().
void
pylith::faults::TestTimeHistorySlipFn::testDbSlipTime(void)
{ // testDbSlipTime
  const char* label = "database ABCD";
  TimeHistorySlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbSlipTime(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbSlipTime);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbSlipTime->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbAmplitude);
  CPPUNIT_ASSERT(0 == slipfn._dbTimeHistory);
} // testDbSlipTime

// ----------------------------------------------------------------------
// Test dbTimeHistory().
void
pylith::faults::TestTimeHistorySlipFn::testDbTimeHistory(void)
{ // testDbTimeHistory
  const char* label = "database ABCDE";
  TimeHistorySlipFn slipfn;
  
  spatialdata::spatialdb::TimeHistory db(label);
  slipfn.dbTimeHistory(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbTimeHistory);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbTimeHistory->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbAmplitude);
  CPPUNIT_ASSERT(0 == slipfn._dbSlipTime);
} // testDbTimeHistory

// ----------------------------------------------------------------------
// Test initialize() in 1-D.
void
pylith::faults::TestTimeHistorySlipFn::testInitialize1D(void)
{ // testInitialize1D
  const char* meshFilename = "data/line2.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/line2_finalslip.spatialdb";
  const char* slipTimeFilename = "data/line2_sliptime.spatialdb";
  const char* timeHistoryFilename = "data/slipfn.timedb";
  const int constraintPts[] = { 3 };
  const double amplitudeE[] = { 2.3 };
  const double slipTimeE[] = { 1.2 };
  const int numConstraintPts = 1;

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
} // testInitialize1D

// ----------------------------------------------------------------------
// Test initialize() in 2-D.
void
pylith::faults::TestTimeHistorySlipFn::testInitialize2D(void)
{ // testInitialize2D
  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* timeHistoryFilename = "data/slipfn.timedb";
  const int constraintPts[] = { 3, 4 };
  const double amplitudeE[] = { 2.3, 0.1, 
				2.4, 0.2};
  const double slipTimeE[] = { 1.2, 1.3 };
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
} // testInitialize2D

// ----------------------------------------------------------------------
// Test initialize() in 3-D.
void
pylith::faults::TestTimeHistorySlipFn::testInitialize3D(void)
{ // testInitialize3D
  const char* meshFilename = "data/tet4.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tet4_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tet4_sliptime.spatialdb";
  const char* timeHistoryFilename = "data/slipfn.timedb";
  const int constraintPts[] = { 3, 4, 5 };
  const double amplitudeE[] = { 2.3, -0.7, 0.1,
				2.4, -0.8, 0.2,
				2.5, -0.9, 0.3 };
  const double slipTimeE[] = { 1.2, 1.3, 1.4 };
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
} // testInitialize3D

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestTimeHistorySlipFn::testSlip(void)
{ // testSlip
  const double slipTimeE[] = { 1.2, 1.3 };
  const double slipE[] = { 0.92, 0.04,
			   0.84, 0.07 };
  const double originTime = 5.064;

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  TimeHistorySlipFn slipfn;
  spatialdata::spatialdb::TimeHistory th;
  _initialize(&mesh, &faultMesh, &slipfn, &th, originTime);
  
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

  const double t = 2.0;
  slipfn.slip(&slip, originTime+t);

  const double tolerance = 1.0e-06;
  int iPoint = 0;
  const ALE::Obj<RealSection>& slipSection = slip.section();
  CPPUNIT_ASSERT(!slipSection.isNull());
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    const int fiberDim = slipSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const double* vals = slipSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vals);
    
    for (int iDim=0; iDim < fiberDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE[iPoint*spaceDim+iDim], 
				   vals[iDim], tolerance);
  } // for
} // testSlip

// ----------------------------------------------------------------------
// Test slipIncr().
void
pylith::faults::TestTimeHistorySlipFn::testSlipIncr(void)
{ // testSlipIncr
  const double slipTimeE[] = { 1.2, 1.3 };
  const double slipE[] = { 0.92, 0.04, 
			   0.984, 0.082};
  const double originTime = 1.064;

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  TimeHistorySlipFn slipfn;
  spatialdata::spatialdb::TimeHistory th;
  _initialize(&mesh, &faultMesh, &slipfn, &th, originTime);

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

  const double t0 = 3.2;
  const double t1 = 9.7;
  slipfn.slipIncr(&slip, originTime+t0, originTime+t1);

  const double tolerance = 1.0e-06;
  int iPoint = 0;
  const ALE::Obj<RealSection>& slipSection = slip.section();
  CPPUNIT_ASSERT(!slipSection.isNull());
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {

    const int fiberDim = slipSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const double* vals = slipSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vals);

    for (int iDim=0; iDim < fiberDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE[iPoint*spaceDim+iDim],
				   vals[iDim], tolerance);
  } // for
} // testSlipIncr

// ----------------------------------------------------------------------
// Initialize TimeHistorySlipFn.
void
pylith::faults::TestTimeHistorySlipFn::_initialize(topology::Mesh* mesh,
					      topology::SubMesh* faultMesh,
					      TimeHistorySlipFn* slipfn,
					      spatialdata::spatialdb::TimeHistory* th,
					      const double originTime)
{ // _initialize
  assert(0 != slipfn);

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
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);

  // Create fault mesh
  int firstFaultVertex = 0;
  int firstLagrangeVertex = mesh->sieveMesh()->getIntSection(faultLabel)->size();
  int firstFaultCell   = mesh->sieveMesh()->getIntSection(faultLabel)->size();
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
} // _initialize

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestTimeHistorySlipFn::_testInitialize(const _TestTimeHistorySlipFn::DataStruct& data)
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
  int firstFaultVertex = 0;
  int firstLagrangeVertex = mesh.sieveMesh()->getIntSection(data.faultLabel)->size();
  int firstFaultCell   = mesh.sieveMesh()->getIntSection(data.faultLabel)->size();
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
  const double originTime = 5.353;
  
  slipfn.initialize(faultMesh, normalizer, originTime);

  const ALE::Obj<SieveMesh::label_sequence>& vertices =
    faultSieveMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

  CPPUNIT_ASSERT(0 != slipfn._parameters);
  const ALE::Obj<RealSection>& finalSlipSection =
    slipfn._parameters->get("slip amplitude").section();
  CPPUNIT_ASSERT(!finalSlipSection.isNull());
  const ALE::Obj<RealSection>& slipTimeSection =
    slipfn._parameters->get("slip time").section();
  CPPUNIT_ASSERT(!slipTimeSection.isNull());

  const double tolerance = 1.0e-06;
  int iPoint = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    CPPUNIT_ASSERT_EQUAL(spaceDim, finalSlipSection->getFiberDimension(*v_iter));
    const double* amplitudeVertex = finalSlipSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != amplitudeVertex);
    for (int iDim=0; iDim < spaceDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.amplitudeE[iPoint*spaceDim+iDim],
				   amplitudeVertex[iDim],
				   tolerance);

    CPPUNIT_ASSERT_EQUAL(1, slipTimeSection->getFiberDimension(*v_iter));
    const double* slipTimeVertex = slipTimeSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != slipTimeVertex);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime,
				 slipTimeVertex[0], tolerance);
  } // for
} // _testInitialize


// End of file 
