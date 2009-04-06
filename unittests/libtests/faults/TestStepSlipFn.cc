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

#include "TestStepSlipFn.hh" // Implementation of class methods

#include "pylith/faults/StepSlipFn.hh" // USES StepSlipFn

#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestStepSlipFn );

// ----------------------------------------------------------------------
namespace pylith {
  namespace faults {
    namespace _TestStepSlipFn {
      struct DataStruct {
	const char* meshFilename;
	const char* faultLabel;
	const int faultId;
	const char* finalSlipFilename;
	const char* slipTimeFilename;
	const int* constraintPts;
	const double* finalSlipE;
	const double* slipTimeE;
	const int numConstraintPts;
      }; // DataStruct
    } // _TestStepSlipFn
  } // faults
} // pylith

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestStepSlipFn::testConstructor(void)
{ // testConstructor
  StepSlipFn slipfn;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbFinalSlip().
void
pylith::faults::TestStepSlipFn::testDbFinalSlip(void)
{ // testDbFinalSlip
  const char* label = "database ABC";
  StepSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbFinalSlip(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbFinalSlip);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbFinalSlip->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbSlipTime);
} // testDbFinalSlip

// ----------------------------------------------------------------------
// Test dbSlipTime().
void
pylith::faults::TestStepSlipFn::testDbSlipTime(void)
{ // testDbSlipTime
  const char* label = "database ABCD";
  StepSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbSlipTime(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbSlipTime);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbSlipTime->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbFinalSlip);
} // testDbSlipTime

// ----------------------------------------------------------------------
// Test initialize() in 1-D.
void
pylith::faults::TestStepSlipFn::testInitialize1D(void)
{ // testInitialize1D
  const char* meshFilename = "data/line2.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/line2_finalslip.spatialdb";
  const char* slipTimeFilename = "data/line2_sliptime.spatialdb";
  const int constraintPts[] = { 3 };
  const double finalSlipE[] = { 2.3 };
  const double slipTimeE[] = { 1.2 };
  const int numConstraintPts = 1;

  _TestStepSlipFn::DataStruct data = {meshFilename,
				      faultLabel,
				      faultId,
				      finalSlipFilename,
				      slipTimeFilename,
				      constraintPts,
				      finalSlipE,
				      slipTimeE,
				      numConstraintPts};
  _testInitialize(data);
} // testInitialize1D

// ----------------------------------------------------------------------
// Test initialize() in 2-D.
void
pylith::faults::TestStepSlipFn::testInitialize2D(void)
{ // testInitialize2D
  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const int constraintPts[] = { 3, 4 };
  const double finalSlipE[] = { 2.3, 0.1, 
				2.4, 0.2 };
  const double slipTimeE[] = { 1.2, 1.3 };
  const int numConstraintPts = 2;

  _TestStepSlipFn::DataStruct data = {meshFilename,
				      faultLabel,
				      faultId,
				      finalSlipFilename,
				      slipTimeFilename,
				      constraintPts,
				      finalSlipE,
				      slipTimeE,
				      numConstraintPts};
  _testInitialize(data);
} // testInitialize2D

// ----------------------------------------------------------------------
// Test initialize() in 3-D.
void
pylith::faults::TestStepSlipFn::testInitialize3D(void)
{ // testInitialize3D
  const char* meshFilename = "data/tet4.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tet4_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tet4_sliptime.spatialdb";
  const int constraintPts[] = { 3, 4, 5 };
  const double finalSlipE[] = { 2.3, -0.7, 0.1,
				2.4, -0.8, 0.2,
				2.5, -0.9, 0.3 };
  const double slipTimeE[] = { 1.2, 1.3, 1.4 };
  const int numConstraintPts = 3;

  _TestStepSlipFn::DataStruct data = {meshFilename,
				      faultLabel,
				      faultId,
				      finalSlipFilename,
				      slipTimeFilename,
				      constraintPts,
				      finalSlipE,
				      slipTimeE,
				      numConstraintPts};
  _testInitialize(data);
} // testInitialize3D

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestStepSlipFn::testSlip(void)
{ // testSlip
  const double slipE[] = { 2.3, 0.1, 
			   0.0, 0.0};
  const double originTime = 5.064;

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  StepSlipFn slipfn;
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

  const double t = 1.234;
  slipfn.slip(slip, originTime+t, faultMesh);

  const double tolerance = 1.0e-06;
  int iPoint = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++iPoint) {
    const int fiberDim = slip->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const real_section_type::value_type* vals = 
      slip->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vals);
    for (int iDim=0; iDim < fiberDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE[iPoint*spaceDim+iDim], vals[iDim], 
				   tolerance);
  } // for
} // testSlip

// ----------------------------------------------------------------------
// Test slipIncr().
void
pylith::faults::TestStepSlipFn::testSlipIncr(void)
{ // testSlipIncr
  const double slipE[] = { 0.0, 0.0, 
			   2.4, 0.2};
  const double originTime = 1.064;

  ALE::Obj<Mesh> faultMesh;
  StepSlipFn slipfn;
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

    const int fiberDim = slip->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const real_section_type::value_type* vals = 
      slip->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vals);

    for (int iDim=0; iDim < fiberDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE[iPoint*spaceDim+iDim], vals[iDim], 
				   tolerance);
  } // for
} // testSlipIncr

// ----------------------------------------------------------------------
// Initialize StepSlipFn.
void
pylith::faults::TestStepSlipFn::_initialize(ALE::Obj<Mesh>* faultMesh,
						 StepSlipFn* slipfn,
						 const double originTime)
{ // _initialize
  assert(0 != slipfn);

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
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
  (*faultMesh)                = new Mesh(mesh->comm(), mesh->getDimension()-1, mesh->debug());
  ALE::Obj<ALE::Mesh> faultBd = NULL;
  CohesiveTopology::createFault(*faultMesh, faultBd,
                                mesh,
                                mesh->getIntSection(faultLabel));
  CohesiveTopology::create(*faultMesh, faultBd, mesh,
                           mesh->getIntSection(faultLabel),
                           faultId,
                           useLagrangeConstraints);
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
  
  spatialdata::units::Nondimensional normalizer;

  // setup StepSlipFn
  slipfn->dbFinalSlip(&dbFinalSlip);
  slipfn->dbSlipTime(&dbSlipTime);
  
  slipfn->initialize(*faultMesh, &cs, normalizer, originTime);
} // _initialize

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestStepSlipFn::_testInitialize(const _TestStepSlipFn::DataStruct& data)
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
  ALE::Obj<Mesh>      faultMesh = new Mesh(mesh->comm(), mesh->getDimension()-1, mesh->debug());
  ALE::Obj<ALE::Mesh> faultBd   = NULL;
  const bool useLagrangeConstraints = true;
  CohesiveTopology::createFault(faultMesh, faultBd,
                                mesh,
                                mesh->getIntSection(data.faultLabel));
  CohesiveTopology::create(faultMesh, faultBd, mesh,
                           mesh->getIntSection(data.faultLabel),
                           data.faultId,
                           useLagrangeConstraints);
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
  
  // setup StepSlipFn
  StepSlipFn slipfn;
  slipfn.dbFinalSlip(&dbFinalSlip);
  slipfn.dbSlipTime(&dbSlipTime);
  
  spatialdata::units::Nondimensional normalizer;
  const double originTime = 5.353;
  
  slipfn.initialize(faultMesh, &cs, normalizer, originTime);

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
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.finalSlipE[iPoint*spaceDim+iDim],
				   vals[iDim],
				   tolerance);

    const double slipTime = vals[spaceDim];

    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime,
				 slipTime, tolerance);
  } // for
} // _testInitialize



// End of file 
