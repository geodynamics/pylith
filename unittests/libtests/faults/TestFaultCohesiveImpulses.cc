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

#include "TestFaultCohesiveImpulses.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveImpulses.hh" // USES FaultCohesiveImpulses

#include "data/CohesiveImpulsesData.hh" // USES CohesiveImpulsesData

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
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveImpulses );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveImpulses::setUp(void)
{ // setUp
  _data = 0;
  _quadrature = new feassemble::Quadrature<topology::SubMesh>();
  CPPUNIT_ASSERT(0 != _quadrature);
  _dbImpulseAmp = 0;
  _flipFault = false;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::faults::TestFaultCohesiveImpulses::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
  delete _dbImpulseAmp; _dbImpulseAmp = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestFaultCohesiveImpulses::testConstructor(void)
{ // testConstructor
  FaultCohesiveImpulses fault;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbImpulseAmp().
void
pylith::faults::TestFaultCohesiveImpulses::testDBImpulseAmp(void)
{ // testDBImpulseAmp
  FaultCohesiveImpulses fault;

  const std::string& label = "test database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  fault.dbImpulseAmp(&db);
  CPPUNIT_ASSERT(fault._dbImpulseAmp);
  CPPUNIT_ASSERT_EQUAL(label, std::string(fault._dbImpulseAmp->label()));
 } // testDBImpulseAmp

// ----------------------------------------------------------------------
// Test threshold().
void
pylith::faults::TestFaultCohesiveImpulses::testThreshold(void)
{ // testThreshold
  FaultCohesiveImpulses fault;

  CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0e-6), fault._threshold); // default

  const PylithScalar value = 1.0e-20;
  fault.threshold(value);
  CPPUNIT_ASSERT_EQUAL(value, fault._threshold);
 } // testThreshold

// ----------------------------------------------------------------------
// Test impulseDOF().
void
pylith::faults::TestFaultCohesiveImpulses::testImpulseDOF(void)
{ // testImpulseDOF
  FaultCohesiveImpulses fault;

  const int ncomps = 2;
  const int dof[ncomps] = { 0, 2 };
  fault.impulseDOF(dof, ncomps);

  CPPUNIT_ASSERT_EQUAL(ncomps, fault.numComponents());
  CPPUNIT_ASSERT_EQUAL(ncomps, int(fault._impulseDOF.size()));
  for (int i=0; i < ncomps; ++i) {
    CPPUNIT_ASSERT_EQUAL(dof[i], fault._impulseDOF[i]);
  } // for
 } // testImpulseDOF

// ----------------------------------------------------------------------
// Test numImpulses().
void
pylith::faults::TestFaultCohesiveImpulses::testNumImpulses(void)
{ // testNumImpulses
  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  FaultCohesiveImpulses fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const int ncomps = 2;
  const int dof[ncomps] = { 0, 2 };
  fault.impulseDOF(dof, ncomps);

  CPPUNIT_ASSERT_EQUAL(_data->numImpulses, fault.numImpulses());
 } // testNumImpulses

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestFaultCohesiveImpulses::testInitialize(void)
{ // testInitialize
  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  FaultCohesiveImpulses fault;
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
    const PylithScalar* orientationVertex =
      orientationSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != orientationVertex);

    const PylithScalar tolerance = 1.0e-06;
    for (int i=0; i < orientationSize; ++i) {
      const int index = iVertex*orientationSize+i;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->orientation[index],
				   orientationVertex[i], tolerance);
    } // for
  } // for

  // Initial tractions
  if (fault._dbImpulseAmp) {
    //fault._fields->get("initial traction").view("INITIAL TRACTIONS"); // DEBUGGING
    const ALE::Obj<RealSection>& amplitudeSection = 
      fault._fields->get("impulse amplitude").section();
    CPPUNIT_ASSERT(!amplitudeSection.isNull());
    const int spaceDim = _data->spaceDim;
    iVertex = 0;
    for (SieveSubMesh::label_sequence::iterator v_iter = verticesBegin;
        v_iter != verticesEnd;
        ++v_iter, ++iVertex) {
      const int fiberDim = amplitudeSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(1, fiberDim);
      const PylithScalar* amplitudeVertex = amplitudeSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(amplitudeVertex);

      const PylithScalar tolerance = 1.0e-06;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->amplitude[iVertex], amplitudeVertex[0], tolerance);
    } // for
  } // if
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::faults::TestFaultCohesiveImpulses::testIntegrateResidual(void)
{ // testIntegrateResidual
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_data->fieldT);
  CPPUNIT_ASSERT(_data->residualIncr);

  topology::Mesh mesh;
  FaultCohesiveImpulses fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const int spaceDim = _data->spaceDim;
  topology::Field<topology::Mesh>& residual = fields.get("residual");
  const ALE::Obj<RealSection>& residualSection = residual.section();
  CPPUNIT_ASSERT(!residualSection.isNull());

  const ALE::Obj<RealSection>& dispSection = fields.get("disp(t)").section();
  CPPUNIT_ASSERT(!dispSection.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  int iVertex = 0;
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex)
    dispSection->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);
  
  const PylithScalar t = 1.0;
  const PylithScalar dt = 1.0;
  fault.timeStep(dt);

  residual.zero();
  { // Integrate residual with disp increment.
    fault.integrateResidual(residual, t, &fields);

    residual.view("RESIDUAL"); // DEBUGGING

    // Check values
    const PylithScalar* valsE = _data->residualIncr;
    iVertex = 0;
    const int fiberDimE = spaceDim;
    const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
    for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
	 v_iter != verticesEnd;
	 ++v_iter, ++iVertex) {
      const int fiberDim = residualSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
      const PylithScalar* vals = residualSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(0 != vals);
      
      for (int i = 0; i < fiberDimE; ++i) {
        const int index = iVertex * spaceDim + i;
        const PylithScalar valE = valsE[index];
        if (fabs(valE) > tolerance)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valE, tolerance);
        else
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
      } // for
    } // for
  } // Integrate residual with disp increment.
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Initialize FaultCohesiveImpulses interface condition.
void
pylith::faults::TestFaultCohesiveImpulses::_initialize(
					topology::Mesh* const mesh,
					FaultCohesiveImpulses* const fault,
					topology::SolutionFields* const fields)
{ // _initialize
  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(fault);
  CPPUNIT_ASSERT(fields);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_quadrature);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);
  
  //mesh->debug(true); // DEBUGGING
  
  _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			  _data->basisDeriv,
			  _data->numQuadPts, _data->numBasis, _data->cellDim,
			  _data->quadPts, _data->numQuadPts, _data->cellDim,
			  _data->quadWts, _data->numQuadPts,
			  _data->spaceDim);
  
  // Setup impulse amplitude
  spatialdata::spatialdb::SimpleDB* db =
      new spatialdata::spatialdb::SimpleDB("impulse amplitude");
  CPPUNIT_ASSERT(db);
  spatialdata::spatialdb::SimpleIOAscii ioImpulseAmp;
  ioImpulseAmp.filename(_data->impulseAmpFilename);
  db->ioHandler(&ioImpulseAmp);
  delete _dbImpulseAmp; _dbImpulseAmp = db;
  fault->dbImpulseAmp(db);

  CPPUNIT_ASSERT(_data->impulseDOF);
  fault->impulseDOF(_data->impulseDOF, _data->numComponents);

  int firstFaultVertex    = 0;
  int firstLagrangeVertex = mesh->sieveMesh()->getIntSection(_data->label)->size();
  int firstFaultCell      = mesh->sieveMesh()->getIntSection(_data->label)->size();
  if (fault->useLagrangeConstraints())
    firstFaultCell += mesh->sieveMesh()->getIntSection(_data->label)->size();
  fault->id(_data->id);
  fault->label(_data->label);
  fault->quadrature(_quadrature);
  
  fault->adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex,
			&firstFaultCell, _flipFault);
  
  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  mesh->nondimensionalize(normalizer);
  
  const PylithScalar upDir[3] = { 0.0, 0.0, 1.0 };
  
  fault->initialize(*mesh, upDir);
  
  // Setup fields
  fields->add("disp(t)", "displacement");
  fields->add("dispIncr(t->t+dt)", "displacement_increment");
  fields->add("velocity(t)", "velocity");
  fields->add("residual", "residual");
  fields->solutionName("dispIncr(t->t+dt)");
  
  const int spaceDim = _data->spaceDim;
  topology::Field<topology::Mesh>& disp = fields->get("disp(t)");
  disp.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  disp.allocate();
  fields->copyLayout("disp(t)");

  fault->verifyConfiguration(*mesh);
} // _initialize

// ----------------------------------------------------------------------
// Determine if vertex is a Lagrange multiplier constraint vertex.
bool
pylith::faults::TestFaultCohesiveImpulses::_isConstraintVertex(const int vertex) const
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
