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

#include "TestStepSlipFn.hh" // Implementation of class methods

#include "pylith/faults/StepSlipFn.hh" // USES StepSlipFn

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
	const PylithScalar* finalSlipE;
	const PylithScalar* slipTimeE;
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
  const PylithScalar finalSlipE[] = { 2.3 };
  const PylithScalar slipTimeE[] = { 1.2 };
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
  const PylithScalar finalSlipE[] = { 2.3, 0.1, 
				2.4, 0.2 };
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
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
  const PylithScalar finalSlipE[] = { 2.3, -0.7, 0.1,
				2.4, -0.8, 0.2,
				2.5, -0.9, 0.3 };
  const PylithScalar slipTimeE[] = { 1.2, 1.3, 1.4 };
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
  const PylithScalar slipE[] = { 2.3, 0.1, 
			   0.0, 0.0};
  const PylithScalar originTime = 5.064;

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  StepSlipFn slipfn;
  _initialize(&mesh, &faultMesh, &slipfn, originTime);
  
  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  CPPUNIT_ASSERT(0 != cs);

  const int spaceDim = cs->spaceDim();
  DM dmMesh = faultMesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  topology::Field<topology::SubMesh> slip(faultMesh);
  slip.newSection(topology::Field<topology::Mesh>::VERTICES_FIELD, spaceDim);
  slip.allocate();

  const PylithScalar t = 1.234;
  slipfn.slip(&slip, originTime+t);

  const PylithScalar tolerance = 1.0e-06;
  PetscScalar       *slipArray;
  int iPoint = 0;
  PetscSection slipSection = slip.petscSection();
  Vec          slipVec     = slip.localVector();
  CPPUNIT_ASSERT(slipSection);CPPUNIT_ASSERT(slipVec);
  err = VecGetArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iPoint) {
    PetscInt dof, off;

    err = PetscSectionGetDof(slipSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipSection, v, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dof);
    for(PetscInt d = 0; d < dof; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE[iPoint*spaceDim+d], slipArray[off+d], tolerance);
  } // for
  err = VecRestoreArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);
} // testSlip

// ----------------------------------------------------------------------
// Initialize StepSlipFn.
void
pylith::faults::TestStepSlipFn::_initialize(topology::Mesh* mesh,
					    topology::SubMesh* faultMesh,
					    StepSlipFn* slipfn,
					    const PylithScalar originTime)
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != faultMesh);
  CPPUNIT_ASSERT(0 != slipfn);
  PetscErrorCode err;

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";

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

  // Create fault mesh
  DM       dmMesh = mesh->dmMesh(), faultBoundaryDM;
  PetscInt firstFaultVertex = 0;
  PetscInt firstLagrangeVertex, firstFaultCell;
  DMLabel  groupField;
  const bool useLagrangeConstraints = true;

  err = DMPlexGetStratumSize(dmMesh, faultLabel, 1, &firstLagrangeVertex);CHECK_PETSC_ERROR(err);
  firstFaultCell = firstLagrangeVertex;
  if (useLagrangeConstraints) {
    firstFaultCell += firstLagrangeVertex;
  }
  err = DMPlexGetLabel(dmMesh, faultLabel, &groupField);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT(groupField);
  ALE::Obj<SieveFlexMesh> faultBoundary = 0;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());CPPUNIT_ASSERT(dmMesh);
  CohesiveTopology::createFault(faultMesh, faultBoundary, faultBoundaryDM,
                                *mesh, groupField);
  CohesiveTopology::create(mesh, *faultMesh, faultBoundary, faultBoundaryDM,
                           groupField,
                           faultId,
                           firstFaultVertex, firstLagrangeVertex, firstFaultCell,
                           useLagrangeConstraints);
  // Need to copy coordinates from mesh to fault mesh since we are not
  // using create() instead of createParallel().
  _setupFaultCoordinates(mesh, faultMesh);

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
  
  slipfn->initialize(*faultMesh, normalizer, originTime);
} // _initialize

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestStepSlipFn::_testInitialize(const _TestStepSlipFn::DataStruct& data)
{ // _testInitialize
  PetscErrorCode err;
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
  DM      dmMesh = mesh.dmMesh(), faultBoundaryDM;
  DMLabel groupField;

  CPPUNIT_ASSERT(!sieveMesh.isNull());CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetLabel(dmMesh, data.faultLabel, &groupField);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT(groupField);
  CohesiveTopology::createFault(&faultMesh, faultBoundary, faultBoundaryDM,
                                mesh, groupField);
  CohesiveTopology::create(&mesh, faultMesh, faultBoundary, faultBoundaryDM,
                           groupField,
                           data.faultId,
                           firstFaultVertex, firstLagrangeVertex, firstFaultCell,
                           useLagrangeConstraints);
  // Need to copy coordinates from mesh to fault mesh since we are
  // using create() instead of createParallel().
  _setupFaultCoordinates(&mesh, &faultMesh);

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
  const PylithScalar originTime = 5.353;
  
  slipfn.initialize(faultMesh, normalizer, originTime);

  CPPUNIT_ASSERT(0 != slipfn._parameters);
  PetscSection finalSlipSection = slipfn._parameters->get("final slip").petscSection();
  Vec          finalSlipVec     = slipfn._parameters->get("final slip").localVector();
  CPPUNIT_ASSERT(finalSlipSection);CPPUNIT_ASSERT(finalSlipVec);
  PetscSection slipTimeSection = slipfn._parameters->get("slip time").petscSection();
  Vec          slipTimeVec     = slipfn._parameters->get("slip time").localVector();
  CPPUNIT_ASSERT(slipTimeSection);CPPUNIT_ASSERT(slipTimeVec);

  const PylithScalar tolerance = 1.0e-06;
  PetscScalar       *finalSlipArray, *slipTimeArray;
  PetscInt           vStart, vEnd, iPoint = 0;
  err = DMPlexGetDepthStratum(faultMesh.dmMesh(), 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = VecGetArray(finalSlipVec, &finalSlipArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(slipTimeVec, &slipTimeArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iPoint) {
    PetscInt dof, off;

    err = PetscSectionGetDof(finalSlipSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(finalSlipSection, v, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dof);
    for(PetscInt d = 0; d < spaceDim; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.finalSlipE[iPoint*spaceDim+d], finalSlipArray[off+d], tolerance);

    err = PetscSectionGetDof(slipTimeSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipTimeSection, v, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(1, dof);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime, slipTimeArray[off], tolerance);
  } // for
  err = VecRestoreArray(finalSlipVec, &finalSlipArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(slipTimeVec, &slipTimeArray);CHECK_PETSC_ERROR(err);
} // _testInitialize

// ----------------------------------------------------------------------
// Setup fault coordinates
void
pylith::faults::TestStepSlipFn::_setupFaultCoordinates(topology::Mesh *mesh, topology::SubMesh *faultMesh)
{ // _setupFaultCoordinates
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh->sieveMesh();
  if (!faultSieveMesh.isNull()) {
    const ALE::Obj<RealSection>& oldCoordSection = mesh->sieveMesh()->getRealSection("coordinates");
    faultSieveMesh->setRealSection("coordinates", oldCoordSection);
  }

  DM              dmMesh      = mesh->dmMesh();
  DM              faultDMMesh = faultMesh->dmMesh();
  const PetscInt  spaceDim    = mesh->dimension();
  IS              subpointIS;
  const PetscInt *points;
  PetscSection    coordSection, fcoordSection;
  PetscInt        vStart, vEnd, ffStart, ffEnd;
  PetscErrorCode  err;

  CPPUNIT_ASSERT(dmMesh);
  CPPUNIT_ASSERT(faultDMMesh);
  err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetHeightStratum(faultDMMesh, 1, &ffStart, &ffEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexCreateSubpointIS(faultDMMesh, &subpointIS);CHECK_PETSC_ERROR(err);
  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMPlexGetCoordinateSection(faultDMMesh, &fcoordSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionSetChart(fcoordSection, vStart, vEnd);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    err = PetscSectionSetDof(fcoordSection, v, spaceDim);CHECK_PETSC_ERROR(err);
  }
  err = PetscSectionSetUp(fcoordSection);CHECK_PETSC_ERROR(err);
  Vec          coordVec, fcoordVec;
  PetscScalar *coords,  *fcoords;
  PetscInt     coordSize;

  err = PetscSectionGetStorageSize(fcoordSection, &coordSize);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  err = VecCreate(mesh->comm(), &fcoordVec);CHECK_PETSC_ERROR(err);
  err = VecSetSizes(fcoordVec, coordSize, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
  err = VecSetFromOptions(fcoordVec);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  err = VecGetArray(fcoordVec, &fcoords);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off, foff;

    // Notice that subpointMap[] does not account for cohesive cells
    err = PetscSectionGetOffset(coordSection, points[v]+(ffEnd-ffStart), &off);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(fcoordSection, v, &foff);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      fcoords[foff+d] = coords[off+d];
    }
    if (!faultSieveMesh.isNull()) {
      const PetscScalar *oldCoords = mesh->sieveMesh()->getRealSection("coordinates")->restrictPoint(points[v]);
      for(PetscInt d = 0; d < spaceDim; ++d) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(oldCoords[d], coords[off+d], 1.0e-6);
      }
    }
  }
  err = ISRestoreIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(fcoordVec, &fcoords);CHECK_PETSC_ERROR(err);
  err = DMSetCoordinatesLocal(faultDMMesh, fcoordVec);CHECK_PETSC_ERROR(err);
  err = VecDestroy(&fcoordVec);CHECK_PETSC_ERROR(err);
} // _setupFaultCoordinates


// End of file 
