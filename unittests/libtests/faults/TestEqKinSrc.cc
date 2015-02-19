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

#include "TestEqKinSrc.hh" // Implementation of class methods

#include "pylith/faults/EqKinSrc.hh" // USES EqKinSrc

#include "TestFaultMesh.hh" // USES createFaultMesh()

#include "pylith/faults/BruneSlipFn.hh" // USES BruneSlipFn
#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestEqKinSrc );

// ----------------------------------------------------------------------
namespace pylith {
  namespace faults {
    namespace _TestEqKinSrc {
      const PylithScalar lengthScale = 1.0e+3;
      const PylithScalar pressureScale = 2.25e+10;
      const PylithScalar timeScale = 1.0;
      const PylithScalar velocityScale = lengthScale / timeScale;
      const PylithScalar densityScale = pressureScale / (velocityScale*velocityScale);
    } // namespace _TestTractPerturbation
  } // faults
} // pylith


// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestEqKinSrc::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  EqKinSrc eqsrc;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test slipFn().
void
pylith::faults::TestEqKinSrc::testSlipFn(void)
{ // testSlipFn
  PYLITH_METHOD_BEGIN;

  BruneSlipFn slipfn;

  EqKinSrc eqsrc;
  eqsrc.slipfn(&slipfn);
  CPPUNIT_ASSERT(&slipfn == eqsrc._slipfn);

  PYLITH_METHOD_END;
} // testSlipFn

// ----------------------------------------------------------------------
// Test initialize(). Use 2-D mesh with Brune slip function to test
// initialize().
void
pylith::faults::TestEqKinSrc::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  topology::Mesh faultMesh;
  EqKinSrc eqsrc;
  BruneSlipFn slipfn;
  const PylithScalar originTime = 2.45;
  _initialize(&mesh, &faultMesh, &eqsrc, &slipfn, originTime);
  
  // Don't have access to details of slip time function, so we can't
  // check parameters. Have to rely on test of slip() for verification
  // of results.

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestEqKinSrc::testSlip(void)
{ // testSlip
  PYLITH_METHOD_BEGIN;

  const PylithScalar finalSlipE[4] = { 2.3, 0.1,   2.4, 0.2,
  };
  const PylithScalar slipTimeE[2] = { 1.2, 1.3 };
  const PylithScalar riseTimeE[2] = { 1.4, 1.5 };
  const PylithScalar originTime = 2.42 / _TestEqKinSrc::timeScale;

  topology::Mesh mesh;
  topology::Mesh faultMesh;
  EqKinSrc eqsrc;
  BruneSlipFn slipfn;
  _initialize(&mesh, &faultMesh, &eqsrc, &slipfn, originTime);
  
  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();CPPUNIT_ASSERT(cs);
  const int spaceDim = cs->spaceDim();

  topology::Field slip(faultMesh);
  slip.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  slip.allocate();

  const PylithScalar t = 2.134 / _TestEqKinSrc::timeScale;
  eqsrc.slip(&slip, originTime+t);

  PetscDM dmMesh = faultMesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  topology::VecVisitorMesh slipVisitor(slip);
  const PetscScalar* slipArray = slipVisitor.localArray();CPPUNIT_ASSERT(slipArray);

  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt v = vStart, iPoint = 0; v < vEnd; ++v, ++iPoint) {
    PylithScalar slipMag = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      slipMag += pow(finalSlipE[iPoint*spaceDim+iDim], 2);
    slipMag = sqrt(slipMag);
    const PylithScalar peakRate = slipMag / riseTimeE[iPoint] * 1.745;
    const PylithScalar tau = slipMag / (exp(1.0) * peakRate);
    const PylithScalar t0 = slipTimeE[iPoint];
    const PylithScalar slipNorm = 1.0 - exp(-(t-t0)/tau) * (1.0 + (t-t0)/tau);

    const PetscInt off = slipVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, slipVisitor.sectionDof(v));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar slipE = finalSlipE[iPoint*spaceDim+d] * slipNorm;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, slipArray[off+d]*_TestEqKinSrc::lengthScale, tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testSlip

// ----------------------------------------------------------------------
// Initialize EqKinSrc.
void
pylith::faults::TestEqKinSrc::_initialize(topology::Mesh* mesh,
					  topology::Mesh* faultMesh,
					  EqKinSrc* eqsrc,
					  BruneSlipFn* slipfn,
					  const PylithScalar originTime)
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(faultMesh);
  CPPUNIT_ASSERT(eqsrc);
  CPPUNIT_ASSERT(slipfn);
  PetscErrorCode err;

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
  const int spaceDim = mesh->dimension();
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);

  // Set scales
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_TestEqKinSrc::lengthScale);
  normalizer.pressureScale(_TestEqKinSrc::pressureScale);
  normalizer.densityScale(_TestEqKinSrc::densityScale);
  normalizer.timeScale(_TestEqKinSrc::timeScale);
  topology::MeshOps::nondimensionalize(mesh, normalizer);

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
  ioRiseTime.filename(peakRateFilename);
  dbRiseTime.ioHandler(&ioRiseTime);

  // setup EqKinSrc
  slipfn->dbFinalSlip(&dbFinalSlip);
  slipfn->dbSlipTime(&dbSlipTime);
  slipfn->dbRiseTime(&dbRiseTime);
  
  eqsrc->originTime(originTime);
  eqsrc->slipfn(slipfn);
  eqsrc->initialize(*faultMesh, normalizer);

  PYLITH_METHOD_END;
} // _initialize


// End of file 
