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

#include "TestConstRateSlipFn.hh" // Implementation of class methods

#include "pylith/faults/ConstRateSlipFn.hh" // USES ConstRateSlipFn

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
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestConstRateSlipFn::testConstructor(void)
{ // testConstructor
  ConstRateSlipFn slipfn;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbSlipRate().
void
pylith::faults::TestConstRateSlipFn::testDbSlipRate(void)
{ // testDbSlipRate
  const char* label = "database ABC";
  ConstRateSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbSlipRate(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbSlipRate);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbSlipRate->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbSlipTime);
} // testDbSlipRate

// ----------------------------------------------------------------------
// Test dbSlipTime().
void
pylith::faults::TestConstRateSlipFn::testDbSlipTime(void)
{ // testDbSlipTime
  const char* label = "database ABCD";
  ConstRateSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbSlipTime(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbSlipTime);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbSlipTime->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbSlipRate);
} // testDbSlipTime

// ----------------------------------------------------------------------
// Test initialize() in 1-D.
void
pylith::faults::TestConstRateSlipFn::testInitialize1D(void)
{ // testInitialize1D
  const char* meshFilename = "data/line2.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* slipRateFilename = "data/line2_sliprate.spatialdb";
  const char* slipTimeFilename = "data/line2_sliptime.spatialdb";
  const int constraintPts[] = { 3 };
  const PylithScalar slipRateE[] = { 0.4 };
  const PylithScalar slipTimeE[] = { 1.2 };
  const int numConstraintPts = 1;

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
} // testInitialize1D

// ----------------------------------------------------------------------
// Test initialize() in 2-D.
void
pylith::faults::TestConstRateSlipFn::testInitialize2D(void)
{ // testInitialize2D
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
} // testInitialize2D

// ----------------------------------------------------------------------
// Test initialize() in 3-D.
void
pylith::faults::TestConstRateSlipFn::testInitialize3D(void)
{ // testInitialize3D
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
} // testInitialize3D

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestConstRateSlipFn::testSlip(void)
{ // testSlip
  const PylithScalar slipRateE[] = { 0.1, 0.2, 
			       0.3, 0.4};
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const PylithScalar originTime = 5.064;

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  ConstRateSlipFn slipfn;
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

  const PylithScalar t = 1.234;
  slipfn.slip(&slip, originTime+t);

  const PylithScalar tolerance = 1.0e-06;
  int iPoint = 0;
  const ALE::Obj<RealSection>& slipSection = slip.section();
  CPPUNIT_ASSERT(!slipSection.isNull());
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    const PylithScalar t0 = slipTimeE[iPoint];
    const int fiberDim = slipSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const PylithScalar* vals = slipSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vals);

    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const PylithScalar slipE = (t - slipTimeE[iPoint]) > 0.0 ?
	slipRateE[iPoint*spaceDim+iDim] * (t - slipTimeE[iPoint]) : 0.0;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, vals[iDim], tolerance);
    } // for
  } // for
} // testSlip

// ----------------------------------------------------------------------
// Initialize ConstRateSlipFn.
void
pylith::faults::TestConstRateSlipFn::_initialize(topology::Mesh* mesh,
						 topology::SubMesh* faultMesh,
					    	 ConstRateSlipFn* slipfn,
						 const PylithScalar originTime)
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != faultMesh);
  CPPUNIT_ASSERT(0 != slipfn);

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
} // _initialize

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestConstRateSlipFn::_testInitialize(const _TestConstRateSlipFn::DataStruct& data)
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

  const ALE::Obj<SieveMesh::label_sequence>& vertices =
    faultSieveMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

  CPPUNIT_ASSERT(0 != slipfn._parameters);
  const ALE::Obj<RealSection>& slipRateSection =
    slipfn._parameters->get("slip rate").section();
  CPPUNIT_ASSERT(!slipRateSection.isNull());
  const ALE::Obj<RealSection>& slipTimeSection =
    slipfn._parameters->get("slip time").section();
  CPPUNIT_ASSERT(!slipTimeSection.isNull());

  const PylithScalar tolerance = 1.0e-06;\
  int iPoint = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    CPPUNIT_ASSERT_EQUAL(spaceDim, slipRateSection->getFiberDimension(*v_iter));
    const PylithScalar* slipRateVertex = slipRateSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != slipRateVertex);
    for (int iDim=0; iDim < spaceDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipRateE[iPoint*spaceDim+iDim],
				   slipRateVertex[iDim],
				   tolerance);

    CPPUNIT_ASSERT_EQUAL(1, slipTimeSection->getFiberDimension(*v_iter));
    const PylithScalar* slipTimeVertex = slipTimeSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != slipTimeVertex);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime,
				 slipTimeVertex[0], tolerance);
  } // for
} // _testInitialize



// End of file 
