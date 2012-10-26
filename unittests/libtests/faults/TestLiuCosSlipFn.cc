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

#include "TestLiuCosSlipFn.hh" // Implementation of class methods

#include "pylith/faults/LiuCosSlipFn.hh" // USES LiuCosSlipFn

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
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestLiuCosSlipFn );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
namespace pylith {
  namespace faults {
    namespace _TestLiuCosSlipFn {
      struct DataStruct {
	const char* meshFilename;
	const char* faultLabel;
	const int faultId;
	const char* finalSlipFilename;
	const char* slipTimeFilename;
	const char* riseTimeFilename;
	const int* constraintPts;
	const PylithScalar* finalSlipE;
	const PylithScalar* slipTimeE;
	const PylithScalar* riseTimeE;
	const int numConstraintPts;
      }; // DataStruct
    } // _TestLiuCosSlipFn
  } // faults
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestLiuCosSlipFn::testConstructor(void)
{ // testConstructor
  LiuCosSlipFn slipfn;
} // testConstructor

// ----------------------------------------------------------------------
// Test dbFinalSlip().
void
pylith::faults::TestLiuCosSlipFn::testDbFinalSlip(void)
{ // testDbFinalSlip
  const char* label = "database ABC";
  LiuCosSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbFinalSlip(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbFinalSlip);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbFinalSlip->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbSlipTime);
  CPPUNIT_ASSERT(0 == slipfn._dbRiseTime);
} // testDbFinalSlip

// ----------------------------------------------------------------------
// Test dbSlipTime().
void
pylith::faults::TestLiuCosSlipFn::testDbSlipTime(void)
{ // testDbSlipTime
  const char* label = "database ABCD";
  LiuCosSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbSlipTime(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbSlipTime);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbSlipTime->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbFinalSlip);
  CPPUNIT_ASSERT(0 == slipfn._dbRiseTime);
} // testDbSlipTime

// ----------------------------------------------------------------------
// Test dbRiseTime().
void
pylith::faults::TestLiuCosSlipFn::testDbRiseTime(void)
{ // testDbRiseTime
  const char* label = "database ABCDE";
  LiuCosSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbRiseTime(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbRiseTime);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbRiseTime->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbFinalSlip);
  CPPUNIT_ASSERT(0 == slipfn._dbSlipTime);
} // testDbRiseTime

// ----------------------------------------------------------------------
// Test initialize() in 1-D.
void
pylith::faults::TestLiuCosSlipFn::testInitialize1D(void)
{ // testInitialize1D
  const char* meshFilename = "data/line2.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/line2_finalslip.spatialdb";
  const char* slipTimeFilename = "data/line2_sliptime.spatialdb";
  const char* riseTimeFilename = "data/line2_risetime.spatialdb";
  const int constraintPts[] = { 3 };
  const PylithScalar finalSlipE[] = { 2.3 };
  const PylithScalar slipTimeE[] = { 1.2 };
  const PylithScalar riseTimeE[] = { 1.4 };
  const int numConstraintPts = 1;

  _TestLiuCosSlipFn::DataStruct data = {meshFilename,
				       faultLabel,
				       faultId,
				       finalSlipFilename,
				       slipTimeFilename,
				       riseTimeFilename,
				       constraintPts,
				       finalSlipE,
				       slipTimeE,
				       riseTimeE,
				       numConstraintPts};
  _testInitialize(data);
} // testInitialize1D

// ----------------------------------------------------------------------
// Test initialize() in 2-D.
void
pylith::faults::TestLiuCosSlipFn::testInitialize2D(void)
{ // testInitialize2D
  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* riseTimeFilename = "data/tri3_risetime.spatialdb";
  const int constraintPts[] = { 3, 4 };
  const PylithScalar finalSlipE[] = { 2.3, 0.1, 
				2.4, 0.2};
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const PylithScalar riseTimeE[] = { 1.4, 1.5 };
  const int numConstraintPts = 2;

  _TestLiuCosSlipFn::DataStruct data = {meshFilename,
				       faultLabel,
				       faultId,
				       finalSlipFilename,
				       slipTimeFilename,
				       riseTimeFilename,
				       constraintPts,
				       finalSlipE,
				       slipTimeE,
				       riseTimeE,
				       numConstraintPts};
  _testInitialize(data);
} // testInitialize2D

// ----------------------------------------------------------------------
// Test initialize() in 3-D.
void
pylith::faults::TestLiuCosSlipFn::testInitialize3D(void)
{ // testInitialize3D
  const char* meshFilename = "data/tet4.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tet4_finalslip.spatialdb";
  const char* slipTimeFilename = "data/tet4_sliptime.spatialdb";
  const char* riseTimeFilename = "data/tet4_risetime.spatialdb";
  const int constraintPts[] = { 3, 4, 5 };
  const PylithScalar finalSlipE[] = { 2.3, -0.7, 0.1,
				2.4, -0.8, 0.2,
				2.5, -0.9, 0.3 };
  const PylithScalar slipTimeE[] = { 1.2, 1.3, 1.4 };
  const PylithScalar riseTimeE[] = { 1.5, 1.6, 1.7 };
  const int numConstraintPts = 3;

  _TestLiuCosSlipFn::DataStruct data = {meshFilename,
				       faultLabel,
				       faultId,
				       finalSlipFilename,
				       slipTimeFilename,
				       riseTimeFilename,
				       constraintPts,
				       finalSlipE,
				       slipTimeE,
				       riseTimeE,
				       numConstraintPts};
  _testInitialize(data);
} // testInitialize3D

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestLiuCosSlipFn::testSlip(void)
{ // testSlip
  const PylithScalar finalSlipE[] = { 2.3, 0.1, 
				0.0, 0.0};
  const PylithScalar slipTimeE[] = { 1.2, 1.3 };
  const PylithScalar riseTimeE[] = { 1.4, 1.5 };
  const PylithScalar originTime = 5.064;

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  LiuCosSlipFn slipfn;
  _initialize(&mesh, &faultMesh, &slipfn, originTime);
  
  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  CPPUNIT_ASSERT(0 != cs);

  const int spaceDim = cs->spaceDim();
  DM dmMesh = faultMesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  topology::Field<topology::SubMesh> slip(faultMesh);
  slip.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  slip.allocate();

  const PylithScalar t = 2.134;
  slipfn.slip(&slip, originTime+t);

  const PylithScalar tolerance = 1.0e-06;
  int iPoint = 0;
  PetscSection slipSection = slip.petscSection();
  Vec          slipVec     = slip.localVector();
  PetscScalar *slipArray;
  CPPUNIT_ASSERT(slipSection);CPPUNIT_ASSERT(slipVec);
  err = VecGetArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iPoint) {
    PylithScalar slipMag = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      slipMag += pow(finalSlipE[iPoint*spaceDim+iDim], 2);
    slipMag = sqrt(slipMag);

    const PylithScalar slipNorm = (slipMag > 0.0) ?
      _slipFn(t - slipTimeE[iPoint], slipMag, riseTimeE[iPoint]) / slipMag :
      0.0;
    PetscInt dof, off;

    err = PetscSectionGetDof(slipSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipSection, v, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dof);

    for(PetscInt d = 0; d < dof; ++d) {
      const PylithScalar slipE = finalSlipE[iPoint*spaceDim+d] * slipNorm;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, slipArray[off+d], tolerance);
    } // for
  } // for
} // testSlip

// ----------------------------------------------------------------------
// Test _slip().
void
pylith::faults::TestLiuCosSlipFn::testSlipTH(void)
{ // testSlipTH
  const PylithScalar t = 0.734;
  const PylithScalar finalSlip = 4.64;
  const PylithScalar riseTime = 3.23;

  const PylithScalar slipE = _slipFn(t, finalSlip, riseTime);

  PylithScalar slip = LiuCosSlipFn::_slipFn(t, finalSlip, riseTime);

  const PylithScalar tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, slip, tolerance);

  slip = LiuCosSlipFn::_slipFn(-0.5, finalSlip, riseTime);
  CPPUNIT_ASSERT_EQUAL(PylithScalar(0.0), slip);

  slip = LiuCosSlipFn::_slipFn(1.0e+10, finalSlip, riseTime);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(finalSlip, slip, tolerance);
} // testSlipTH

// ----------------------------------------------------------------------
// Initialize LiuCosSlipFn.
void
pylith::faults::TestLiuCosSlipFn::_initialize(topology::Mesh* mesh,
					      topology::SubMesh* faultMesh,
					      LiuCosSlipFn* slipfn,
					      const PylithScalar originTime)
{ // _initialize
  assert(0 != slipfn);

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* finalSlipFilename = "data/tri3_finalslipB.spatialdb";
  const char* slipTimeFilename = "data/tri3_sliptime.spatialdb";
  const char* riseTimeFilename = "data/tri3_risetime.spatialdb";

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
  int firstFaultVertex    = 0;
  int firstLagrangeVertex = mesh->sieveMesh()->getIntSection(faultLabel)->size();
  int firstFaultCell      = mesh->sieveMesh()->getIntSection(faultLabel)->size();
  const bool useLagrangeConstraints = true;
  if (useLagrangeConstraints) {
    firstFaultCell += mesh->sieveMesh()->getIntSection(faultLabel)->size();
  }
  ALE::Obj<SieveFlexMesh> faultBoundary = 0;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CohesiveTopology::createFault(faultMesh, faultBoundary,
                                *mesh, sieveMesh->getIntSection(faultLabel));
  CohesiveTopology::create(mesh, *faultMesh, faultBoundary, 
                           sieveMesh->getIntSection(faultLabel),
                           faultId, firstFaultVertex, firstLagrangeVertex,
			   firstFaultCell, useLagrangeConstraints);
  // Need to copy coordinates from mesh to fault mesh since we are not
  // using create() instead of createParallel().
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh->sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  const ALE::Obj<RealSection>& oldCoordSection = sieveMesh->getRealSection("coordinates");
  faultSieveMesh->setRealSection("coordinates", oldCoordSection);
  DM              dmMesh = faultMesh->dmMesh();
  IS              subpointMap;
  const PetscInt *points;
  PetscSection    coordSection;
  PetscInt        vStart, vEnd;
  PetscErrorCode  err;
  CPPUNIT_ASSERT(dmMesh);

  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMComplexGetSubpointMap(dmMesh, &subpointMap);CHECK_PETSC_ERROR(err);
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
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
  err = ISGetIndices(subpointMap, &points);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(coordSection, v, &off);CHECK_PETSC_ERROR(err);
    const PetscScalar *oldCoords = oldCoordSection->restrictPoint(points[v]);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      coords[off+d] = oldCoords[d];
    }
  }
  err = ISRestoreIndices(subpointMap, &points);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  err = DMSetCoordinatesLocal(dmMesh, coordVec);CHECK_PETSC_ERROR(err);

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
  ioRiseTime.filename(riseTimeFilename);
  dbRiseTime.ioHandler(&ioRiseTime);

  spatialdata::units::Nondimensional normalizer;

  // setup LiuCosSlipFn
  slipfn->dbFinalSlip(&dbFinalSlip);
  slipfn->dbSlipTime(&dbSlipTime);
  slipfn->dbRiseTime(&dbRiseTime);
  
  slipfn->initialize(*faultMesh, normalizer, originTime);
} // _initialize

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestLiuCosSlipFn::_testInitialize(const _TestLiuCosSlipFn::DataStruct& data)
{ // _testInitialize
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
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CohesiveTopology::createFault(&faultMesh, faultBoundary,
                                mesh, sieveMesh->getIntSection(data.faultLabel));
  CohesiveTopology::create(&mesh, faultMesh, faultBoundary, 
                           sieveMesh->getIntSection(data.faultLabel),
                           data.faultId, firstFaultVertex, firstLagrangeVertex,
			   firstFaultCell, useLagrangeConstraints);
  // Need to copy coordinates from mesh to fault mesh since we are not
  // using create() instead of createParallel().
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh.sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  const ALE::Obj<RealSection>& oldCoordSection = sieveMesh->getRealSection("coordinates");
  faultSieveMesh->setRealSection("coordinates", oldCoordSection);
  DM              dmMesh = faultMesh.dmMesh();
  IS              subpointMap;
  const PetscInt *points;
  PetscSection    coordSection;
  PetscInt        vStart, vEnd;
  PetscErrorCode  err;
  CPPUNIT_ASSERT(dmMesh);

  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMComplexGetSubpointMap(dmMesh, &subpointMap);CHECK_PETSC_ERROR(err);
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
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
  err = ISGetIndices(subpointMap, &points);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(coordSection, v, &off);CHECK_PETSC_ERROR(err);
    const PetscScalar *oldCoords = oldCoordSection->restrictPoint(points[v]);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      coords[off+d] = oldCoords[d];
    }
  }
  err = ISRestoreIndices(subpointMap, &points);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  err = DMSetCoordinatesLocal(dmMesh, coordVec);CHECK_PETSC_ERROR(err);

  // Setup databases
  spatialdata::spatialdb::SimpleDB dbFinalSlip("final slip");
  spatialdata::spatialdb::SimpleIOAscii ioFinalSlip;
  ioFinalSlip.filename(data.finalSlipFilename);
  dbFinalSlip.ioHandler(&ioFinalSlip);
  
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
  ioSlipTime.filename(data.slipTimeFilename);
  dbSlipTime.ioHandler(&ioSlipTime);
  
  spatialdata::spatialdb::SimpleDB dbRiseTime("rise time");
  spatialdata::spatialdb::SimpleIOAscii ioRiseTime;
  ioRiseTime.filename(data.riseTimeFilename);
  dbRiseTime.ioHandler(&ioRiseTime);

  // setup LiuCosSlipFn
  LiuCosSlipFn slipfn;
  slipfn.dbFinalSlip(&dbFinalSlip);
  slipfn.dbSlipTime(&dbSlipTime);
  slipfn.dbRiseTime(&dbRiseTime);
  
  spatialdata::units::Nondimensional normalizer;
  const PylithScalar originTime = 5.353;
  
  slipfn.initialize(faultMesh, normalizer, originTime);

  CPPUNIT_ASSERT(0 != slipfn._parameters);
  PetscSection finalSlipSection = slipfn._parameters->get("final slip").petscSection();
  Vec          finalSlipVec     = slipfn._parameters->get("final slip").localVector();
  PetscScalar *finalSlipArray;
  CPPUNIT_ASSERT(finalSlipSection);CPPUNIT_ASSERT(finalSlipVec);
  PetscSection slipTimeSection = slipfn._parameters->get("slip time").petscSection();
  Vec          slipTimeVec     = slipfn._parameters->get("slip time").localVector();
  PetscScalar *slipTimeArray;
  CPPUNIT_ASSERT(slipTimeSection);CPPUNIT_ASSERT(slipTimeVec);
  PetscSection riseTimeSection = slipfn._parameters->get("rise time").petscSection();
  Vec          riseTimeVec     = slipfn._parameters->get("rise time").localVector();
  PetscScalar *riseTimeArray;
  CPPUNIT_ASSERT(riseTimeSection);CPPUNIT_ASSERT(riseTimeVec);

  const PylithScalar tolerance = 1.0e-06;
  int iPoint = 0;
  err = VecGetArray(finalSlipVec, &finalSlipArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(slipTimeVec,  &slipTimeArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(riseTimeVec,  &riseTimeArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iPoint) {
    PetscInt fsdof, fsoff, stdof, stoff, rtdof, rtoff;

    err = PetscSectionGetDof(finalSlipSection, v, &fsdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(finalSlipSection, v, &fsoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(slipTimeSection, v, &stdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipTimeSection, v, &stoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(riseTimeSection, v, &rtdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(riseTimeSection, v, &rtoff);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fsdof);
    for(PetscInt d = 0; d < fsdof; ++d)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(data.finalSlipE[iPoint*spaceDim+d], finalSlipArray[fsoff+d], tolerance);

    CPPUNIT_ASSERT_EQUAL(1, stdof);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.slipTimeE[iPoint]+originTime, slipTimeArray[stoff], tolerance);

    CPPUNIT_ASSERT_EQUAL(1, rtdof);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(data.riseTimeE[iPoint], riseTimeArray[rtoff], tolerance);
  } // for
  err = VecRestoreArray(finalSlipVec, &finalSlipArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(slipTimeVec,  &slipTimeArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(riseTimeVec,  &riseTimeArray);CHECK_PETSC_ERROR(err);
} // _testInitialize

// ----------------------------------------------------------------------
// Slip time function.
PylithScalar
pylith::faults::TestLiuCosSlipFn::_slipFn(const PylithScalar t,
					  const PylithScalar finalSlip,
					  const PylithScalar riseTime)
{ // _slipFn
  const float tau = riseTime * 1.525;
  const float tau1 = 0.13 * tau;
  const float tau2 = tau - tau1;
  const float Cn = 
    M_PI /  (1.4 * M_PI * tau1 + 1.2 * tau1 + 0.3 * M_PI * tau2);
  
  PylithScalar slip = 0.0;
  if (t <= tau1) {
    slip = 0.7*t - 0.7*tau1/M_PI*sin(M_PI*t/tau1) 
      - 0.6*tau1/(0.5*M_PI)*(cos(0.5*M_PI*t/tau1) - 1.0);
    slip *= Cn;
  } else if (t <= 2.0*tau1) {
    slip = 1.0*t - 0.7*tau1/M_PI*sin(M_PI*t/tau1) + 
      0.3*tau2/M_PI*sin(M_PI*(t-tau1)/tau2) + 
      1.2*tau1/M_PI - 0.3*tau1;
    slip *= Cn;
  } else if (t <= tau) {
    slip = 0.3*t + 0.3*tau2/M_PI*sin(M_PI*(t-tau1)/tau2) + 
      1.1*tau1 + 1.2*tau1/M_PI;
    slip *= Cn;	
  } else
    slip = 1.0;
  slip *= finalSlip;

  return slip;
} // _slipFn


// End of file 
