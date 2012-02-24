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

#include "TestEqKinSrc.hh" // Implementation of class methods

#include "pylith/faults/EqKinSrc.hh" // USES EqKinSrc


#include "pylith/faults/BruneSlipFn.hh" // USES BruneSlipFn
#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestEqKinSrc );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestEqKinSrc::testConstructor(void)
{ // testConstructor
  EqKinSrc eqsrc;
} // testConstructor

// ----------------------------------------------------------------------
// Test slipFn().
void
pylith::faults::TestEqKinSrc::testSlipFn(void)
{ // testSlipFn
  BruneSlipFn slipfn;

  EqKinSrc eqsrc;
  eqsrc.slipfn(&slipfn);
  CPPUNIT_ASSERT(&slipfn == eqsrc._slipfn);
} // testSlipFn

// ----------------------------------------------------------------------
// Test initialize(). Use 2-D mesh with Brune slip function to test
// initialize().
void
pylith::faults::TestEqKinSrc::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  EqKinSrc eqsrc;
  BruneSlipFn slipfn;
  const PylithScalar originTime = 2.45;
  _initialize(&mesh, &faultMesh, &eqsrc, &slipfn, originTime);
  
  // Don't have access to details of slip time function, so we can't
  // check parameters. Have to rely on test of slip() for verification
  // of results.
} // testInitialize

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestEqKinSrc::testSlip(void)
{ // testSlip
  const PylithScalar finalSlipE[] = { 2.3, 0.1, 
				2.4, 0.2};
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const PylithScalar riseTimeE[] = { 1.4, 1.5 };
  const PylithScalar originTime = 2.42;

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  EqKinSrc eqsrc;
  BruneSlipFn slipfn;
  _initialize(&mesh, &faultMesh, &eqsrc, &slipfn, originTime);
  
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
  eqsrc.slip(&slip, originTime+t);

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
    const PylithScalar tau = slipMag / (exp(1.0) * peakRate);
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
pylith::faults::TestEqKinSrc::testSlipIncr(void)
{ // testSlip
  const PylithScalar finalSlipE[] = { 2.3, 0.1, 
				2.4, 0.2};
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const PylithScalar riseTimeE[] = { 1.4, 1.5 };
  const PylithScalar originTime = -4.29;

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  EqKinSrc eqsrc;
  BruneSlipFn slipfn;
  _initialize(&mesh, &faultMesh, &eqsrc, &slipfn, originTime);
  
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
  const PylithScalar t1 = 2.525;
  eqsrc.slipIncr(&slip, originTime+t0, originTime+t1);

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
    const PylithScalar tau = slipMag / (exp(1.0) * peakRate);
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
	finalSlipE[iPoint*spaceDim+iDim] * (slipNorm1 - slipNorm0);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, vals[iDim], tolerance);
    } // for
  } // for
} // testSlipIncr

// ----------------------------------------------------------------------
// Initialize EqKinSrc.
void
pylith::faults::TestEqKinSrc::_initialize(topology::Mesh* mesh,
					  topology::SubMesh* faultMesh,
					  EqKinSrc* eqsrc,
					  BruneSlipFn* slipfn,
					  const PylithScalar originTime)
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != faultMesh);
  CPPUNIT_ASSERT(0 != eqsrc);
  CPPUNIT_ASSERT(0 != slipfn);

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* peakRateFilename = "data/tri3_risetime.spatialdb";

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
  ioRiseTime.filename(peakRateFilename);
  dbRiseTime.ioHandler(&ioRiseTime);

  spatialdata::units::Nondimensional normalizer;

  // setup EqKinSrc
  slipfn->dbFinalSlip(&dbFinalSlip);
  slipfn->dbSlipTime(&dbSlipTime);
  slipfn->dbRiseTime(&dbRiseTime);
  
  eqsrc->originTime(originTime);
  eqsrc->slipfn(slipfn);
  eqsrc->initialize(*faultMesh, normalizer);
} // _initialize


// End of file 
