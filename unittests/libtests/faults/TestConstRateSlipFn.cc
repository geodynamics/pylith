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

#include "TestConstRateSlipFn.hh" // Implementation of class methods

#include "pylith/faults/ConstRateSlipFn.hh" // USES ConstRateSlipFn

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
	const PylithScalar* slipRateE;
	const PylithScalar* slipTimeE;
	const int numConstraintPts;
      }; // DataStruct
    } // _TestConstRateSlipFn
  } // faults
} // pylith

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestConstRateSlipFn::testConstructor(void)
{ // testConstructor
  ConstRateSlipFn slipfn;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbSlipRate().
void
pylith::faults::TestConstRateSlipFn::testDbSlipRate(void)
{ // testDbSlipRate
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
  const PylithScalar slipRateE[] = { 0.4 };
  const PylithScalar slipTimeE[] = { 1.2 };
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
  const PylithScalar slipRateE[] = { 0.1, 0.2, 
			       0.3, 0.4 };
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
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
  const PylithScalar slipRateE[] = { 1.6, -0.7, 0.1,
			       1.7, -0.8, 0.2,
			       1.8, -0.9, 0.3 };
  const PylithScalar slipTimeE[] = { 1.2, 1.3, 1.4 };
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
  const PylithScalar slipRateE[] = { 0.1, 0.2, 
			       0.3, 0.4};
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const PylithScalar originTime = 5.064;

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  ConstRateSlipFn slipfn;
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
  int iPoint = 0;
  PetscSection slipSection = slip.petscSection();
  Vec          slipVec     = slip.localVector();
  PetscScalar *slipArray;
  CPPUNIT_ASSERT(slipSection);CPPUNIT_ASSERT(slipVec);
  err = VecGetArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iPoint) {
    const PylithScalar t0 = slipTimeE[iPoint];
    PetscInt dof, off;

    err = PetscSectionGetDof(slipSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipSection, v, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dof);
    for(PetscInt d = 0; d < dof; ++d) {
      const PylithScalar slipE = (t - t0) > 0.0 ? slipRateE[iPoint*spaceDim+d] * (t - t0) : 0.0;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, slipArray[off+d], tolerance);
    }
  } // for
  err = VecRestoreArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);
} // testSlip

// ----------------------------------------------------------------------
// Initialize ConstRateSlipFn.
void
pylith::faults::TestConstRateSlipFn::_initialize(topology::Mesh* mesh,
						 topology::SubMesh* faultMesh,
					    	 ConstRateSlipFn* slipfn,
						 const PylithScalar originTime)
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != faultMesh);
  CPPUNIT_ASSERT(0 != slipfn);
  PetscErrorCode err;

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* slipRateFilename = "data/tri3_sliprate.spatialdb";
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
  DM       dmMesh = mesh->dmMesh();
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
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CohesiveTopology::createFault(faultMesh, faultBoundary,
                                *mesh, groupField);
  CohesiveTopology::create(mesh, *faultMesh, faultBoundary, 
                           sieveMesh->getIntSection(faultLabel),
                           faultId,
                           firstFaultVertex, firstLagrangeVertex, firstFaultCell,
                           useLagrangeConstraints);
  // Need to copy coordinates from mesh to fault mesh since we are not
  // using create() instead of createParallel().
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh->sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  const ALE::Obj<RealSection>& oldCoordSection = sieveMesh->getRealSection("coordinates");
  faultSieveMesh->setRealSection("coordinates", oldCoordSection);
  DM              faultDMMesh = faultMesh->dmMesh();
  IS              subpointIS;
  const PetscInt *points;
  PetscSection    coordSection;
  PetscInt        vStart, vEnd;

  CPPUNIT_ASSERT(faultDMMesh);
  err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexCreateSubpointIS(faultDMMesh, &subpointIS);CHECK_PETSC_ERROR(err);
  err = DMPlexGetCoordinateSection(faultDMMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionSetChart(coordSection, vStart, vEnd);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    err = PetscSectionSetDof(coordSection, v, spaceDim);CHECK_PETSC_ERROR(err);
  }
  err = PetscSectionSetUp(coordSection);CHECK_PETSC_ERROR(err);
  Vec          coordVec;
  PetscScalar *coords;
  PetscInt     coordSize;

  err = PetscSectionGetStorageSize(coordSection, &coordSize);CHECK_PETSC_ERROR(err);
  err = VecCreate(mesh->comm(), &coordVec);CHECK_PETSC_ERROR(err);
  err = VecSetSizes(coordVec, coordSize, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
  err = VecSetFromOptions(coordVec);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(coordSection, v, &off);CHECK_PETSC_ERROR(err);
    const PetscScalar *oldCoords = oldCoordSection->restrictPoint(points[v]);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      coords[off+d] = oldCoords[d];
    }
  }
  err = ISRestoreIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  err = DMSetCoordinatesLocal(faultDMMesh, coordVec);CHECK_PETSC_ERROR(err);

  // Setup databases
  spatialdata::spatialdb::SimpleDB dbSlipRate("slip rate");
  spatialdata::spatialdb::SimpleIOAscii ioSlipRate;
  ioSlipRate.filename(slipRateFilename);
  dbSlipRate.ioHandler(&ioSlipRate);
  
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
  ioSlipTime.filename(slipTimeFilename);
  dbSlipTime.ioHandler(&ioSlipTime);
  
  spatialdata::units::Nondimensional normalizer;

  // setup ConstRateSlipFn
  slipfn->dbSlipRate(&dbSlipRate);
  slipfn->dbSlipTime(&dbSlipTime);
  
  slipfn->initialize(*faultMesh, normalizer, originTime);
} // _initialize

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestConstRateSlipFn::_testInitialize(const _TestConstRateSlipFn::DataStruct& data)
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
  DM       dmMesh = mesh.dmMesh();
  PetscInt firstFaultVertex = 0;
  PetscInt firstLagrangeVertex, firstFaultCell;
  DMLabel  groupField;
  const bool useLagrangeConstraints = true;

  err = DMPlexGetStratumSize(dmMesh, data.faultLabel, 1, &firstLagrangeVertex);CHECK_PETSC_ERROR(err);
  firstFaultCell = firstLagrangeVertex;
  if (useLagrangeConstraints) {
    firstFaultCell += firstLagrangeVertex;
  }
  err = DMPlexGetLabel(dmMesh, data.faultLabel, &groupField);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT(groupField);
  ALE::Obj<SieveFlexMesh> faultBoundary = 0;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CohesiveTopology::createFault(&faultMesh, faultBoundary,
                                mesh, groupField);
  CohesiveTopology::create(&mesh, faultMesh, faultBoundary, 
                           sieveMesh->getIntSection(data.faultLabel),
                           data.faultId,
                           firstFaultVertex, firstLagrangeVertex, firstFaultCell,
                           useLagrangeConstraints);
  // Need to copy coordinates from mesh to fault mesh since we are not
  // using create() instead of createParallel().
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh.sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  const ALE::Obj<RealSection>& oldCoordSection = sieveMesh->getRealSection("coordinates");
  faultSieveMesh->setRealSection("coordinates", oldCoordSection);
  DM              faultDMMesh = faultMesh.dmMesh();
  IS              subpointIS;
  const PetscInt *points;
  PetscSection    coordSection;
  PetscInt        vStart, vEnd;

  CPPUNIT_ASSERT(faultDMMesh);
  err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexCreateSubpointIS(faultDMMesh, &subpointIS);CHECK_PETSC_ERROR(err);
  err = DMPlexGetCoordinateSection(faultDMMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionSetChart(coordSection, vStart, vEnd);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    err = PetscSectionSetDof(coordSection, v, spaceDim);CHECK_PETSC_ERROR(err);
  }
  err = PetscSectionSetUp(coordSection);CHECK_PETSC_ERROR(err);
  Vec          coordVec;
  PetscScalar *coords;
  PetscInt     coordSize;

  err = PetscSectionGetStorageSize(coordSection, &coordSize);CHECK_PETSC_ERROR(err);
  err = VecCreate(mesh.comm(), &coordVec);CHECK_PETSC_ERROR(err);
  err = VecSetSizes(coordVec, coordSize, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
  err = VecSetFromOptions(coordVec);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(coordSection, v, &off);CHECK_PETSC_ERROR(err);
    const PetscScalar *oldCoords = oldCoordSection->restrictPoint(points[v]);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      coords[off+d] = oldCoords[d];
    }
  }
  err = ISRestoreIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  err = DMSetCoordinatesLocal(faultDMMesh, coordVec);CHECK_PETSC_ERROR(err);

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
  
  spatialdata::units::Nondimensional normalizer;
  const PylithScalar originTime = 5.353;
  
  slipfn.initialize(faultMesh, normalizer, originTime);

  CPPUNIT_ASSERT(0 != slipfn._parameters);
  PetscSection slipRateSection = slipfn._parameters->get("slip rate").petscSection();
  Vec          slipRateVec     = slipfn._parameters->get("slip rate").localVector();
  PetscScalar *slipRateArray;
  CPPUNIT_ASSERT(slipRateSection);CPPUNIT_ASSERT(slipRateVec);
  PetscSection slipTimeSection = slipfn._parameters->get("slip time").petscSection();
  Vec          slipTimeVec     = slipfn._parameters->get("slip time").localVector();
  PetscScalar *slipTimeArray;
  CPPUNIT_ASSERT(slipTimeSection);CPPUNIT_ASSERT(slipTimeVec);

  const PylithScalar tolerance = 1.0e-06;
  int iPoint = 0;
  err = VecGetArray(slipRateVec, &slipRateArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(slipTimeVec, &slipTimeArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iPoint) {
    PetscInt srdof, sroff, stdof, stoff;

    err = PetscSectionGetDof(slipRateSection, v, &srdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipRateSection, v, &sroff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(slipTimeSection, v, &stdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipTimeSection, v, &stoff);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, srdof);
    for(PetscInt d = 0; d < spaceDim; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipRateE[iPoint*spaceDim+d], slipRateArray[sroff+d], tolerance);

    CPPUNIT_ASSERT_EQUAL(1, stdof);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime, slipTimeArray[stoff], tolerance);
  } // for
  err = VecRestoreArray(slipRateVec, &slipRateArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(slipTimeVec, &slipTimeArray);CHECK_PETSC_ERROR(err);
} // _testInitialize



// End of file 
