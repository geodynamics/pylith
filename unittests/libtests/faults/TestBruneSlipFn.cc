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
  const double finalSlipE[] = { 2.3, 0.1, 
				0.0, 0.0};
  const double slipTimeE[] = { 1.2, 1.3 };
  const double peakRateE[] = { 1.4, 1.5 };
  const double originTime = 5.064;

  ALE::Obj<Mesh> faultMesh;
  BruneSlipFn slipfn;
  _initialize(&faultMesh, &slipfn, originTime);
  
  const int spaceDim = faultMesh->getDimension() + 1;

  const ALE::Obj<Mesh::label_sequence>& vertices = faultMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();
  ALE::Obj<real_section_type> slip = 
    new real_section_type(faultMesh->comm(), faultMesh->debug());
  slip->setChart(real_section_type::chart_type(*std::min_element(vertices->begin(), 
								 vertices->end()), 
					       *std::max_element(vertices->begin(), vertices->end())+1));
  slip->setFiberDimension(vertices, spaceDim);
  faultMesh->allocate(slip);
  CPPUNIT_ASSERT(!slip.isNull());

  const double t = 2.134;
  slipfn.slip(slip, originTime+t, faultMesh);

  const double tolerance = 1.0e-06;
  int iPoint = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    double slipMag = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      slipMag += pow(finalSlipE[iPoint*spaceDim+iDim], 2);
    slipMag = sqrt(slipMag);
    const double tau = 
      (slipMag > 0.0) ? slipMag / (exp(1.0) * peakRateE[iPoint]) : 1.0;
    const double t0 = slipTimeE[iPoint];
    const double slipNorm = 1.0 - exp(-(t-t0)/tau) * (1.0 + (t-t0)/tau);
    const int fiberDim = slip->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const real_section_type::value_type* vals = 
      slip->restrictPoint(*v_iter);
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
pylith::faults::TestBruneSlipFn::testSlipIncr(void)
{ // testSlipIncr
  const double finalSlipE[] = { 2.3, 0.1, 
				0.0, 0.0};
  const double slipTimeE[] = { 1.2, 1.3 };
  const double peakRateE[] = { 1.4, 1.5 };
  const double originTime = 1.064;

  ALE::Obj<Mesh> faultMesh;
  BruneSlipFn slipfn;
  _initialize(&faultMesh, &slipfn, originTime);

  const int spaceDim = faultMesh->getDimension() + 1;

  const ALE::Obj<Mesh::label_sequence>& vertices = faultMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();
  ALE::Obj<real_section_type> slip =
    new real_section_type(faultMesh->comm(), faultMesh->debug());
  slip->setChart(real_section_type::chart_type(*std::min_element(vertices->begin(), 
								 vertices->end()), 
					       *std::max_element(vertices->begin(), vertices->end())+1));
  slip->setFiberDimension(vertices, spaceDim);
  faultMesh->allocate(slip);
  CPPUNIT_ASSERT(!slip.isNull());

  const double t0 = 1.234;
  const double t1 = 3.635;
  slipfn.slipIncr(slip, originTime+t0, originTime+t1, faultMesh);

  const double tolerance = 1.0e-06;
  int iPoint = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    double slipMag = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      slipMag += pow(finalSlipE[iPoint*spaceDim+iDim], 2);
    slipMag = sqrt(slipMag);
    const double tau = 
      (slipMag > 0.0) ? slipMag / (exp(1.0) * peakRateE[iPoint]) : 1.0;
    const double tRef = slipTimeE[iPoint];
    const double slipNorm0 = 
      (t0 > tRef) ? 1.0 - exp(-(t0-tRef)/tau) * (1.0 + (t0-tRef)/tau) : 0.0;
    const double slipNorm1 =
      (t1 > tRef) ? 1.0 - exp(-(t1-tRef)/tau) * (1.0 + (t1-tRef)/tau) : 0.0;

    const int fiberDim = slip->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const real_section_type::value_type* vals = 
      slip->restrictPoint(*v_iter);
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
pylith::faults::TestBruneSlipFn::testSlipTH(void)
{ // testSlipTH
  const double t = 0.734;
  const double finalSlip = 4.64;
  const double peakRate = 3.23;

  const double tau = finalSlip / (exp(1.0) * peakRate);
  const double slipE = finalSlip * (1.0 - exp(-t/tau) * (1.0 + t/tau));

  double slip = BruneSlipFn::_slipFn(t, finalSlip, peakRate);

  const double tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, slip, tolerance);

  slip = BruneSlipFn::_slipFn(-0.5, finalSlip, peakRate);
  CPPUNIT_ASSERT_EQUAL(0.0, slip);

  slip = BruneSlipFn::_slipFn(1.0e+10, finalSlip, peakRate);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(finalSlip, slip, tolerance);
} // testSlipTH

// ----------------------------------------------------------------------
// Initialize BruneSlipFn.
void
pylith::faults::TestBruneSlipFn::_initialize(ALE::Obj<Mesh>* faultMesh,
					     BruneSlipFn* slipfn,
					     const double originTime)
{ // _initialize
  assert(0 != slipfn);

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslipB.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* peakRateFilename = "data/tri3_peakrate.spatialdb";

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
  const bool useLagrangeConstraints = true;
  CohesiveTopology::create(faultMesh, mesh, 
			   mesh->getIntSection(faultLabel),
			   faultId);
  CPPUNIT_ASSERT(!faultMesh->isNull());
  // Need to copy coordinates from mesh to fault mesh since we are not
  // using create() instead of createParallel().
  (*faultMesh)->setRealSection("coordinates", 
			       mesh->getRealSection("coordinates"));

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
  slipfn->dbFinalSlip(&dbFinalSlip);
  slipfn->dbSlipTime(&dbSlipTime);
  slipfn->dbPeakRate(&dbPeakRate);
  
  slipfn->initialize(*faultMesh, &cs, originTime);
} // _initialize

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
  // Need to copy coordinates from mesh to fault mesh since we are not
  // using create() instead of createParallel().
  faultMesh->setRealSection("coordinates", 
			    mesh->getRealSection("coordinates"));

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
  BruneSlipFn slipfn;
  slipfn.dbFinalSlip(&dbFinalSlip);
  slipfn.dbSlipTime(&dbSlipTime);
  slipfn.dbPeakRate(&dbPeakRate);
  
  const double originTime = 5.353;
  
  slipfn.initialize(faultMesh, &cs, originTime);

  const double tolerance = 1.0e-06;

  const ALE::Obj<Mesh::label_sequence>& vertices = 
    faultMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  int iPoint = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    const int fiberDim = slipfn._parameters->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim+2, fiberDim);
    
    const real_section_type::value_type* vals = 
      slipfn._parameters->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vals);

    for (int iDim=0; iDim < spaceDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.finalSlipE[iPoint*spaceDim+iDim],
				   vals[iDim],
				   tolerance);

    const double peakRate = vals[spaceDim  ];
    const double slipTime = vals[spaceDim+1];

    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.peakRateE[iPoint], peakRate, tolerance);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime,
				 slipTime, tolerance);
  } // for
} // _testInitialize



// End of file 
