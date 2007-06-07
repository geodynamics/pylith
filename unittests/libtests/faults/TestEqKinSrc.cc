// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestEqKinSrc.hh" // Implementation of class methods

#include "pylith/faults/EqKinSrc.hh" // USES EqKinSrc

#include "pylith/faults/BruneSlipFn.hh" // USES BruneSlipFn
#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestEqKinSrc );

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
  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* peakRateFilename = "data/tri3_peakrate.spatialdb";
  const int constraintPts[] = { 3, 4 };
  const double finalSlipE[] = { 2.3, 0.1, 
				2.4, 0.2};
  const double slipTimeE[] = { 1.2, 1.3 };
  const double peakRateE[] = { 1.4, 1.5 };
  const int numConstraintPts = 2;

  typedef std::set<Mesh::point_type>::const_iterator vert_iterator;  

  // Setup mesh
  ALE::Obj<Mesh> mesh;
  meshio::MeshIOAscii meshIO;
  meshIO.filename(meshFilename);
  meshIO.debug(false);
  meshIO.interpolate(false);
  meshIO.read(&mesh);
  CPPUNIT_ASSERT(!mesh.isNull());
  const int spaceDim = mesh->getDimension();
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);

  // Create fault mesh
  ALE::Obj<Mesh> faultMesh;
  const bool useLagrangeConstraints = true;
  CohesiveTopology::create(&faultMesh, mesh, 
			   mesh->getIntSection(faultLabel),
			   faultId);
  CPPUNIT_ASSERT(!faultMesh.isNull());

  // Create set of constraint vertices
  std::set<Mesh::point_type> eqsrcVertices;
  for (int i=0; i < numConstraintPts; ++i)
    eqsrcVertices.insert(constraintPts[i]);
  CPPUNIT_ASSERT_EQUAL(numConstraintPts, int(eqsrcVertices.size()));
  
  // Setup databases
  spatialdata::spatialdb::SimpleDB dbFinalSlip("final slip");
  spatialdata::spatialdb::SimpleIOAscii ioFinalSlip;
  ioFinalSlip.filename(finalSlipFilename);
  dbFinalSlip.ioHandler(&ioFinalSlip);
  
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
  ioSlipTime.filename(slipTimeFilename);
  dbSlipTime.ioHandler(&ioSlipTime);
  
  spatialdata::spatialdb::SimpleDB dbPeakRate("peak rate");
  spatialdata::spatialdb::SimpleIOAscii ioPeakRate;
  ioPeakRate.filename(peakRateFilename);
  dbPeakRate.ioHandler(&ioPeakRate);

  // setup EqKinSrc
  BruneSlipFn slipFn;
  slipFn.dbFinalSlip(&dbFinalSlip);
  slipFn.dbSlipTime(&dbSlipTime);
  slipFn.dbPeakRate(&dbPeakRate);
  
  EqKinSrc eqsrc;
  eqsrc.slipfn(&slipFn);
  eqsrc.initialize(mesh, faultMesh, eqsrcVertices, &cs);

  // Don't have access to details of slip time function, so we can't
  // check parameters. Have to rely on test of slip() for verification
  // of results.
} // testInitialize

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestEqKinSrc::testSlip(void)
{ // testSlip
  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* peakRateFilename = "data/tri3_peakrate.spatialdb";
  const int constraintPts[] = { 3, 4 };
  const double finalSlipE[] = { 2.3, 0.1, 
				2.4, 0.2};
  const double slipTimeE[] = { 1.2, 1.3 };
  const double peakRateE[] = { 1.4, 1.5 };
  const int numConstraintPts = 2;

  ALE::Obj<Mesh> mesh;
  meshio::MeshIOAscii meshIO;
  meshIO.filename(meshFilename);
  meshIO.debug(false);
  meshIO.interpolate(false);
  meshIO.read(&mesh);
  CPPUNIT_ASSERT(!mesh.isNull());
  const int spaceDim = mesh->getDimension();
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);

  // Create fault mesh
  ALE::Obj<Mesh> faultMesh;
  const bool useLagrangeConstraints = true;
  CohesiveTopology::create(&faultMesh, mesh, 
			   mesh->getIntSection(faultLabel),
			   faultId);
  CPPUNIT_ASSERT(!faultMesh.isNull());

  // Create set of constraint vertices
  std::set<Mesh::point_type> eqsrcVertices;
  for (int i=0; i < numConstraintPts; ++i)
    eqsrcVertices.insert(constraintPts[i]);
  CPPUNIT_ASSERT_EQUAL(numConstraintPts, int(eqsrcVertices.size()));
  
  // Setup databases
  spatialdata::spatialdb::SimpleDB dbFinalSlip("final slip");
  spatialdata::spatialdb::SimpleIOAscii ioFinalSlip;
  ioFinalSlip.filename(finalSlipFilename);
  dbFinalSlip.ioHandler(&ioFinalSlip);
  
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
  ioSlipTime.filename(slipTimeFilename);
  dbSlipTime.ioHandler(&ioSlipTime);
  
  spatialdata::spatialdb::SimpleDB dbPeakRate("peak rate");
  spatialdata::spatialdb::SimpleIOAscii ioPeakRate;
  ioPeakRate.filename(peakRateFilename);
  dbPeakRate.ioHandler(&ioPeakRate);

  // setup EqKinSrc
  BruneSlipFn slipFn;
  slipFn.dbFinalSlip(&dbFinalSlip);
  slipFn.dbSlipTime(&dbSlipTime);
  slipFn.dbPeakRate(&dbPeakRate);
  
  EqKinSrc eqsrc;
  eqsrc.slipfn(&slipFn);
  eqsrc.initialize(mesh, faultMesh, eqsrcVertices, &cs);
  
  const double t = 2.134;
  const ALE::Obj<real_section_type>& slip = eqsrc.slip(t, eqsrcVertices);
  CPPUNIT_ASSERT(!slip.isNull());

  int iPoint = 0;
  const double tolerance = 1.0e-06;
  typedef std::set<Mesh::point_type>::const_iterator vert_iterator;  
  const vert_iterator vBegin = eqsrcVertices.begin();
  const vert_iterator vEnd = eqsrcVertices.end();
  for (vert_iterator v_iter=vBegin; v_iter != vEnd; ++v_iter, ++iPoint) {
    double slipMag = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      slipMag += pow(finalSlipE[iPoint*spaceDim+iDim], 2);
    slipMag = sqrt(slipMag);
    const double tau = slipMag / (exp(1.0) * peakRateE[iPoint]);
    const double t0 = slipTimeE[iPoint];
    const double slipNorm = 1.0 - exp(-(t-t0)/tau) * (1.0 + (t-t0)/tau);

    const int fiberDim = slip->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const real_section_type::value_type* vals = 
      slip->restrictPoint(*v_iter);
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const double slipE = finalSlipE[iPoint*spaceDim+iDim] * slipNorm;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, vals[iDim], tolerance);
    } // for
  } // for
} // testSlip


// End of file 
