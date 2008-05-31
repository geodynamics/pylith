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

#include "TestConstRateSlipFn.hh" // Implementation of class methods

#include "pylith/faults/ConstRateSlipFn.hh" // USES ConstRateSlipFn

#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

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
	const double* slipRateE;
	const double* slipTimeE;
	const int numConstraintPts;
      }; // DataStruct
    } // _TestConstRateSlipFn
  } // faults
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestConstRateSlipFn::testConstructor(void)
{ // testConstructor
  ConstRateSlipFn slipfn;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbFinalSlip().
void
pylith::faults::TestConstRateSlipFn::testDbSlipRate(void)
{ // testDbFinalSlip
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
  const double slipRateE[] = { 0.4 };
  const double slipTimeE[] = { 1.2 };
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
  const double slipRateE[] = { 0.1, 0.2, 
			       0.3, 0.4 };
  const double slipTimeE[] = { 1.2, 1.3 };
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
  const double slipRateE[] = { 1.6, -0.7, 0.1,
			       1.7, -0.8, 0.2,
			       1.8, -0.9, 0.3 };
  const double slipTimeE[] = { 1.2, 1.3, 1.4 };
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
  const double slipRateE[] = { 0.1, 0.2, 
			       0.3, 0.4};
  const double slipTimeE[] = { 1.2, 1.3 };
  const double originTime = 5.064;

  ALE::Obj<Mesh> faultMesh;
  ConstRateSlipFn slipfn;
  _initialize(&faultMesh, &slipfn, originTime);
  
  const int spaceDim = faultMesh->getDimension() + 1;

  const double t = 2.134;
  const ALE::Obj<real_section_type>& slip = slipfn.slip(originTime+t, faultMesh);
  CPPUNIT_ASSERT(!slip.isNull());

  const double tolerance = 1.0e-06;
  
  const ALE::Obj<Mesh::label_sequence>& vertices = 
    faultMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  int iPoint = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    const double t0 = slipTimeE[iPoint];
    const int fiberDim = slip->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const real_section_type::value_type* vals = 
      slip->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vals);

    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const double slipE = 
	slipRateE[iPoint*spaceDim+iDim] * (t - slipTimeE[iPoint]);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, vals[iDim], tolerance);
    } // for
  } // for
} // testSlip

// ----------------------------------------------------------------------
// Test slipIncr().
void
pylith::faults::TestConstRateSlipFn::testSlipIncr(void)
{ // testSlipIncr
  const double slipRateE[] = { 0.1, 0.2, 
			       0.3, 0.4 };
  const double slipTimeE[] = { 1.2, 1.3 };
  const double originTime = 1.064;

  ALE::Obj<Mesh> faultMesh;
  ConstRateSlipFn slipfn;
  _initialize(&faultMesh, &slipfn, originTime);

  const int spaceDim = faultMesh->getDimension() + 1;

  const double t0 = 1.234;
  const double t1 = 3.635;
  const ALE::Obj<real_section_type>& slip = 
    slipfn.slipIncr(originTime+t0, originTime+t1, faultMesh);
  CPPUNIT_ASSERT(!slip.isNull());

  const double tolerance = 1.0e-06;

  const ALE::Obj<Mesh::label_sequence>& vertices = 
    faultMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  int iPoint = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {

    const int fiberDim = slip->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const real_section_type::value_type* vals = 
      slip->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vals);

    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const double tRef = (slipTimeE[iPoint] > t0) ? slipTimeE[iPoint] : t0;
      const double slipE = 
	slipRateE[iPoint*spaceDim+iDim] * (t1 - tRef);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, vals[iDim], tolerance);
    } // for
  } // for
} // testSlipIncr

// ----------------------------------------------------------------------
// Initialize ConstRateSlipFn.
void
pylith::faults::TestConstRateSlipFn::_initialize(ALE::Obj<Mesh>* faultMesh,
						 ConstRateSlipFn* slipfn,
						 const double originTime)
{ // _initialize
  assert(0 != slipfn);

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* slipRateFilename = "data/tri3_sliprate.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";

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
  spatialdata::spatialdb::SimpleDB dbSlipRate("slip rate");
  spatialdata::spatialdb::SimpleIOAscii ioSlipRate;
  ioSlipRate.filename(slipRateFilename);
  dbSlipRate.ioHandler(&ioSlipRate);
  
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
  ioSlipTime.filename(slipTimeFilename);
  dbSlipTime.ioHandler(&ioSlipTime);
  
  // setup ConstRateSlipFn
  slipfn->dbSlipRate(&dbSlipRate);
  slipfn->dbSlipTime(&dbSlipTime);
  
  slipfn->initialize(*faultMesh, &cs, originTime);
} // _initialize

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestConstRateSlipFn::_testInitialize(const _TestConstRateSlipFn::DataStruct& data)
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
    CPPUNIT_ASSERT_EQUAL(spaceDim+1, fiberDim);
    
    const real_section_type::value_type* vals = 
      slipfn._parameters->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vals);

    for (int iDim=0; iDim < spaceDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipRateE[iPoint*spaceDim+iDim],
				   vals[iDim],
				   tolerance);

    const double slipTime = vals[spaceDim];

    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime,
				 slipTime, tolerance);
  } // for
} // _testInitialize



// End of file 
