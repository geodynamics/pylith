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

#include "TestBruneSlipFn.hh" // Implementation of class methods

#include "pylith/faults/BruneSlipFn.hh" // USES BruneSlipFn

#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

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
	const char* peakRateFilename;
	const int* constraintPts;
	const double* finalSlipE;
	const double* slipTimeE;
	const double* peakRateE;
	const int numConstraintPts;
      }; // DataStruct
    } // _TestBruneSlipFn
  } // faults
} // pylith

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
  CPPUNIT_ASSERT(0 == slipfn._dbPeakRate);
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
  CPPUNIT_ASSERT(0 == slipfn._dbPeakRate);
} // testDbSlipTime

// ----------------------------------------------------------------------
// Test dbPeakRate().
void
pylith::faults::TestBruneSlipFn::testDbPeakRate(void)
{ // testDbPeakRate
  const char* label = "database ABCDE";
  BruneSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbPeakRate(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbPeakRate);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbPeakRate->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbFinalSlip);
  CPPUNIT_ASSERT(0 == slipfn._dbSlipTime);
} // testDbPeakRate

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
  const char* peakRateFilename = "data/line2_peakrate.spatialdb";
  const int constraintPts[] = { 3 };
  const double finalSlipE[] = { 2.3 };
  const double slipTimeE[] = { 1.2 };
  const double peakRateE[] = { 1.4 };
  const int numConstraintPts = 1;

  _TestBruneSlipFn::DataStruct data = {meshFilename,
				       faultLabel,
				       faultId,
				       finalSlipFilename,
				       slipTimeFilename,
				       peakRateFilename,
				       constraintPts,
				       finalSlipE,
				       slipTimeE,
				       peakRateE,
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
  const char* peakRateFilename = "data/tri3_peakrate.spatialdb";
  const int constraintPts[] = { 3, 4 };
  const double finalSlipE[] = { 2.3, 0.1, 
				2.4, 0.2};
  const double slipTimeE[] = { 1.2, 1.3 };
  const double peakRateE[] = { 1.4, 1.5 };
  const int numConstraintPts = 2;

  _TestBruneSlipFn::DataStruct data = {meshFilename,
				       faultLabel,
				       faultId,
				       finalSlipFilename,
				       slipTimeFilename,
				       peakRateFilename,
				       constraintPts,
				       finalSlipE,
				       slipTimeE,
				       peakRateE,
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
  const char* peakRateFilename = "data/tet4_peakrate.spatialdb";
  const int constraintPts[] = { 3, 4, 5 };
  const double finalSlipE[] = { 2.3, -0.7, 0.1,
				2.4, -0.8, 0.2,
				2.5, -0.9, 0.3 };
  const double slipTimeE[] = { 1.2, 1.3, 1.4 };
  const double peakRateE[] = { 1.5, 1.6, 1.7 };
  const int numConstraintPts = 3;

  _TestBruneSlipFn::DataStruct data = {meshFilename,
				       faultLabel,
				       faultId,
				       finalSlipFilename,
				       slipTimeFilename,
				       peakRateFilename,
				       constraintPts,
				       finalSlipE,
				       slipTimeE,
				       peakRateE,
				       numConstraintPts};
  _testInitialize(data);
} // testInitialize3D

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestBruneSlipFn::testSlip(void)
{ // testSlip
  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslipB.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* peakRateFilename = "data/tri3_peakrate.spatialdb";
  const int constraintPts[] = { 3, 4 };
  const double finalSlipE[] = { 2.3, 0.1, 
				0.0, 0.0};
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

  // setup BruneSlipFn
  BruneSlipFn slipFn;
  slipFn.dbFinalSlip(&dbFinalSlip);
  slipFn.dbSlipTime(&dbSlipTime);
  slipFn.dbPeakRate(&dbPeakRate);
  
  slipFn.initialize(mesh, faultMesh, eqsrcVertices, &cs);
  
  const double t = 2.134;
  const ALE::Obj<real_section_type>& slip = slipFn.slip(t, eqsrcVertices);
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

// ----------------------------------------------------------------------
// Test slipIncr().
void
pylith::faults::TestBruneSlipFn::testSlipIncr(void)
{ // testSlipIncr
  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslipB.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* peakRateFilename = "data/tri3_peakrate.spatialdb";
  const int constraintPts[] = { 3, 4 };
  const double finalSlipE[] = { 2.3, 0.1, 
				0.0, 0.0};
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

  // setup BruneSlipFn
  BruneSlipFn slipFn;
  slipFn.dbFinalSlip(&dbFinalSlip);
  slipFn.dbSlipTime(&dbSlipTime);
  slipFn.dbPeakRate(&dbPeakRate);
  
  slipFn.initialize(mesh, faultMesh, eqsrcVertices, &cs);
  
  const double t0 = 1.234;
  const double t1 = 3.635;
  const ALE::Obj<real_section_type>& slip = 
    slipFn.slipIncr(t0, t1, eqsrcVertices);
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
    const double tRef = slipTimeE[iPoint];
    const double slipNorm0 = 
      (t0 > tRef) ? 1.0 - exp(-(t0-tRef)/tau) * (1.0 + (t0-tRef)/tau) : 0.0;
    const double slipNorm1 =
      (t1 > tRef) ? 1.0 - exp(-(t1-tRef)/tau) * (1.0 + (t1-tRef)/tau) : 0.0;

    const int fiberDim = slip->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const real_section_type::value_type* vals = 
      slip->restrictPoint(*v_iter);
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
pylith::faults::TestBruneSlipFn::testSlipTH(void)
{ // testSlipTH
  const double t = 0.734;
  const double finalSlip = 4.64;
  const double peakRate = 3.23;

  const double tau = finalSlip / (exp(1.0) * peakRate);
  const double slipE = finalSlip * (1.0 - exp(-t/tau) * (1.0 + t/tau));

  double slip = BruneSlipFn::_slip(t, finalSlip, peakRate);

  const double tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, slip, tolerance);

  slip = BruneSlipFn::_slip(-0.5, finalSlip, peakRate);
  CPPUNIT_ASSERT_EQUAL(0.0, slip);

  slip = BruneSlipFn::_slip(1.0e+10, finalSlip, peakRate);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(finalSlip, slip, tolerance);
} // testSlipTH

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestBruneSlipFn::_testInitialize(const _TestBruneSlipFn::DataStruct& data)
{ // _testInitialize
  typedef std::set<Mesh::point_type>::const_iterator vert_iterator;  

  // Setup mesh
  ALE::Obj<Mesh> mesh;
  meshio::MeshIOAscii meshIO;
  meshIO.filename(data.meshFilename);
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
			   mesh->getIntSection(data.faultLabel),
			   data.faultId);
  CPPUNIT_ASSERT(!faultMesh.isNull());

  // Create set of constraint vertices
  std::set<Mesh::point_type> eqsrcVertices;
  for (int i=0; i < data.numConstraintPts; ++i)
    eqsrcVertices.insert(data.constraintPts[i]);
  CPPUNIT_ASSERT_EQUAL(data.numConstraintPts, int(eqsrcVertices.size()));
  
  // Setup databases
  spatialdata::spatialdb::SimpleDB dbFinalSlip("final slip");
  spatialdata::spatialdb::SimpleIOAscii ioFinalSlip;
  ioFinalSlip.filename(data.finalSlipFilename);
  dbFinalSlip.ioHandler(&ioFinalSlip);
  
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
  ioSlipTime.filename(data.slipTimeFilename);
  dbSlipTime.ioHandler(&ioSlipTime);
  
  spatialdata::spatialdb::SimpleDB dbPeakRate("peak rate");
  spatialdata::spatialdb::SimpleIOAscii ioPeakRate;
  ioPeakRate.filename(data.peakRateFilename);
  dbPeakRate.ioHandler(&ioPeakRate);

  // setup BruneSlipFn
  BruneSlipFn slipFn;
  slipFn.dbFinalSlip(&dbFinalSlip);
  slipFn.dbSlipTime(&dbSlipTime);
  slipFn.dbPeakRate(&dbPeakRate);
  
  slipFn.initialize(mesh, faultMesh, eqsrcVertices, &cs);

  // Check parameter sections
  const double tolerance = 1.0e-06;
  CPPUNIT_ASSERT(0 != slipFn._parameters);
  const ALE::Obj<real_section_type>& finalSlip = 
    slipFn._parameters->getReal("final slip");
  const ALE::Obj<real_section_type>& slipTime = 
    slipFn._parameters->getReal("slip time");
  const ALE::Obj<real_section_type>& peakRate = 
    slipFn._parameters->getReal("peak rate");
  CPPUNIT_ASSERT(!finalSlip.isNull());
  CPPUNIT_ASSERT(!slipTime.isNull());
  CPPUNIT_ASSERT(!peakRate.isNull());

  int iPoint = 0;
  const vert_iterator vBegin = eqsrcVertices.begin();
  const vert_iterator vEnd = eqsrcVertices.end();
  for (vert_iterator v_iter=vBegin; v_iter != vEnd; ++v_iter, ++iPoint) {
    { // final slip
      const int fiberDim = finalSlip->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
      const real_section_type::value_type* vals = 
	finalSlip->restrictPoint(*v_iter);
      for (int iDim=0; iDim < fiberDim; ++iDim)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.finalSlipE[iPoint*spaceDim+iDim],
				     vals[iDim],
				     tolerance);
    } // final slip
    
    { // slip time
      const int fiberDim = slipTime->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(1, fiberDim);
      const real_section_type::value_type* vals = 
	slipTime->restrictPoint(*v_iter);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint], vals[0], tolerance);
    } // slip time

    { // peak rate
      const int fiberDim = peakRate->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(1, fiberDim);
      const real_section_type::value_type* vals = 
	peakRate->restrictPoint(*v_iter);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.peakRateE[iPoint], vals[0], tolerance);
    } // peak rate
  } // for
} // _testInitialize



// End of file 
