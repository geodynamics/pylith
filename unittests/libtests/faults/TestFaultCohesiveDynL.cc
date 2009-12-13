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

#include "TestFaultCohesiveDynL.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveDynL.hh" // USES FaultCohesiveDynL

#include "data/CohesiveDynLData.hh" // USES CohesiveDynLData

#include "pylith/faults/EqKinSrc.hh" // USES EqKinSrc
#include "pylith/faults/BruneSlipFn.hh" // USES BruneSlipFn
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

#include <stdexcept> // USES runtime_error

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveDynL );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveDynL::setUp(void)
{ // setUp
  _data = 0;
  _quadrature = new feassemble::Quadrature<topology::SubMesh>();
  CPPUNIT_ASSERT(0 != _quadrature);
  _dbInitialTract = 0;

  _flipFault = false;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::faults::TestFaultCohesiveDynL::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
  delete _dbInitialTract; _dbInitialTract = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestFaultCohesiveDynL::testConstructor(void)
{ // testConstructor
  FaultCohesiveDynL fault;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbInitialTract().
void
pylith::faults::TestFaultCohesiveDynL::testDBInitialTract(void)
{ // testDBInitialTract
  FaultCohesiveDynL fault;

  const std::string& label = "test database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  fault.dbInitialTract(&db);
  CPPUNIT_ASSERT(0 != fault._dbInitialTract);
  CPPUNIT_ASSERT_EQUAL(label, std::string(fault._dbInitialTract->label()));
 } // testDBInitialTract

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestFaultCohesiveDynL::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  FaultCohesiveDynL fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const ALE::Obj<SieveSubMesh>& faultSieveMesh = fault._faultMesh->sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  SieveSubMesh::renumbering_type& renumbering = 
    faultSieveMesh->getRenumbering();
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = 
    faultSieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  int iVertex = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    CPPUNIT_ASSERT(renumbering.find(_data->constraintVertices[iVertex]) !=
		   renumbering.end());
    CPPUNIT_ASSERT_EQUAL(renumbering[_data->constraintVertices[iVertex]],
			 *v_iter);
  } // for
  CPPUNIT_ASSERT_EQUAL(_data->numConstraintVert, iVertex);

  // Check orientation
  //fault._fields->get("orientation").view("ORIENTATION"); // DEBUGGING
  const ALE::Obj<RealSection>& orientationSection = 
    fault._fields->get("orientation").section();
  CPPUNIT_ASSERT(!orientationSection.isNull());
  const int spaceDim = _data->spaceDim;
  const int orientationSize = spaceDim*spaceDim;
  iVertex = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    const int fiberDim = orientationSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(orientationSize, fiberDim);
    const double* orientationVertex =
      orientationSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != orientationVertex);

    const double tolerance = 1.0e-06;
    for (int i=0; i < orientationSize; ++i) {
      const int index = iVertex*orientationSize+i;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->orientation[index],
				   orientationVertex[i], tolerance);
    } // for
  } // for

  // Check area
  const ALE::Obj<RealSection>& areaSection =
    fault._fields->get("area").section();
  CPPUNIT_ASSERT(!areaSection.isNull());
  iVertex = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    const int fiberDim = areaSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(1, fiberDim);
    const double* areaVertex = areaSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != areaVertex);

    const double tolerance = 1.0e-06;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->area[iVertex], areaVertex[0],
				 tolerance);
  } // for
} // testInitialize

// ----------------------------------------------------------------------
// Test constrainSolnSpace() for sticking case.
void
pylith::faults::TestFaultCohesiveDynL::testConstrainSolnSpaceStick(void)
{ // testConstrainSolnSpaceStick
  // STUFF GOES HERE
} // testConstrainSolnSpaceStick

// ----------------------------------------------------------------------
// Test constrainSolnSpace() for slipping case.
void
pylith::faults::TestFaultCohesiveDynL::testConstrainSolnSpaceSlip(void)
{ // testConstrainSolnSpaceSlip
  // STUFF GOES HERE
} // testConstrainSolnSpaceSlip

// ----------------------------------------------------------------------
// Test constrainSolnSpace() for opening case.
void
pylith::faults::TestFaultCohesiveDynL::testConstrainSolnSpaceOpen(void)
{ // testConstrainSolnSpaceOpen
  // STUFF GOES HERE
} // testConstrainSolnSpaceOpen

// ----------------------------------------------------------------------
// Test updateStateVars().
void
pylith::faults::TestFaultCohesiveDynL::testUpdateStateVars(void)
{ // testUpdateStateVars
  // STUFF GOES HERE
} // testUpdateStateVars

// ----------------------------------------------------------------------
// Test calcTractions().
void
pylith::faults::TestFaultCohesiveDynL::testCalcTractions(void)
{ // testCalcTractions
  topology::Mesh mesh;
  FaultCohesiveDynL fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  
  // STUFF GOES HERE
} // testCalcTractions

// ----------------------------------------------------------------------
// Initialize FaultCohesiveDynL interface condition.
void
pylith::faults::TestFaultCohesiveDynL::_initialize(
					topology::Mesh* const mesh,
					FaultCohesiveDynL* const fault,
					topology::SolutionFields* const fields)
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != fault);
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != _quadrature);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);
  
  //mesh->debug(true); // DEBUGGING
  
  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  mesh->nondimensionalize(normalizer);
  
  _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			  _data->basisDeriv,
			  _data->numQuadPts, _data->numBasis, _data->cellDim,
			  _data->quadPts, _data->numQuadPts, _data->cellDim,
			  _data->quadWts, _data->numQuadPts,
			  _data->spaceDim);
  
  // Setup initial tractions
  spatialdata::spatialdb::SimpleDB* db =
      new spatialdata::spatialdb::SimpleDB("initial tractions");
  CPPUNIT_ASSERT(0 != db);
  spatialdata::spatialdb::SimpleIOAscii ioInitialTract;
  ioInitialTract.filename(_data->initialTractFilename);
  db->ioHandler(&ioInitialTract);
  delete _dbInitialTract; _dbInitialTract = db;

  int firstFaultVertex    = 0;
  int firstLagrangeVertex = mesh->sieveMesh()->getIntSection(_data->label)->size();
  int firstFaultCell      = mesh->sieveMesh()->getIntSection(_data->label)->size();
  if (fault->useLagrangeConstraints()) {
    firstFaultCell += mesh->sieveMesh()->getIntSection(_data->label)->size();
  }
  fault->id(_data->id);
  fault->label(_data->label);
  fault->quadrature(_quadrature);
  
  fault->adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell, _flipFault);
  
  const double upDir[] = { 0.0, 0.0, 1.0 };
  const double normalDir[] = { 1.0, 0.0, 0.0 };
  
  fault->initialize(*mesh, upDir, normalDir); 
  
  // Setup fields
  fields->add("residual", "residual");
  fields->add("disp(t)", "displacement");
  fields->add("dispIncr(t->t+dt)", "displacement_increment");
  fields->solutionName("dispIncr(t->t+dt)");
  
  const int spaceDim = _data->spaceDim;
  topology::Field<topology::Mesh>& residual = fields->get("residual");
  residual.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  residual.allocate();
  fields->copyLayout("residual");
} // _initialize

// ----------------------------------------------------------------------
// Determine if vertex is a Lagrange multiplier constraint vertex.
bool
pylith::faults::TestFaultCohesiveDynL::_isConstraintVertex(const int vertex) const
{ // _isConstraintVertex
  assert(0 != _data);

  const int numConstraintVert = _data->numConstraintVert;
  bool isFound = false;
  for (int i=0; i < _data->numConstraintVert; ++i)
    if (_data->constraintVertices[i] == vertex) {
      isFound = true;
      break;
    } // if
  return isFound;
} // _isConstraintVertex


// End of file 
