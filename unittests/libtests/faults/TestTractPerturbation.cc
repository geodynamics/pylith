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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestTractPerturbation.hh" // Implementation of class methods

#include "pylith/faults/TractPerturbation.hh" // USES TractPerturbation

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
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestTractPerturbation );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestTractPerturbation::testConstructor(void)
{ // testConstructor
  TractPerturbation eqsrc;
} // testConstructor

// ----------------------------------------------------------------------
// Test label().
void
pylith::faults::TestTractPerturbation::testLabel(void)
{ // testLabel
  const std::string& label = "nucleation";

  TractPerturbation tract;
  tract.label(label.c_str());
  CPPUNIT_ASSERT_EQUAL(label, tract._label);
} // testLabel

// ----------------------------------------------------------------------
// Test hasParameter().
void
pylith::faults::TestTractPerturbation::testHasParameter(void)
{ // testHasParameter
  spatialdata::spatialdb::SimpleDB db;
  
  TractPerturbation tract;

  // no values
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_initial_value"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_rate_of_change"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_change_in_value"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_rate_start_time"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_change_start_time"));

  // initial value
  tract.dbInitial(&db);
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_initial_value"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_rate_of_change"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_change_in_value"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_rate_start_time"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_change_start_time"));

  // change value
  tract.dbChange(&db);
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_initial_value"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_rate_of_change"));
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_change_in_value"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_rate_start_time"));
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_change_start_time"));

  // rate value, remove change
  tract.dbRate(&db);
  tract.dbChange(0);
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_initial_value"));
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_rate_of_change"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_change_in_value"));
  CPPUNIT_ASSERT_EQUAL(true, tract.hasParameter("traction_rate_start_time"));
  CPPUNIT_ASSERT_EQUAL(false, tract.hasParameter("traction_change_start_time"));
} // testHasParameter

// ----------------------------------------------------------------------
// Test initialize() using 2-D mesh.
void
pylith::faults::TestTractPerturbation::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  TractPerturbation tract;
  _initialize(&mesh, &faultMesh, &tract);
  
  // Rely on testTraction() for verification of results.
} // testInitialize

// ----------------------------------------------------------------------
// Test calculate() using 2-D mesh().
void
pylith::faults::TestTractPerturbation::testCalculate(void)
{ // testCalculate
  const PylithScalar tractionE[4] = { 
    -1.0*(-2.0+1.0), -1.0*(1.0-0.5), // initial + change
    -1.0*(-2.1), -1.0*(1.1), // initial
  };

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  TractPerturbation tract;
  _initialize(&mesh, &faultMesh, &tract);
  
  const PylithScalar t = 2.134;
  tract.calculate(t);


  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  CPPUNIT_ASSERT(cs);
  const int spaceDim = cs->spaceDim();

  DM             faultDMMesh = faultMesh.dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  CPPUNIT_ASSERT(faultDMMesh);
  err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);

  //traction.view("TRACTION"); // DEBUGGING

  const PylithScalar tolerance = 1.0e-06;
  int iPoint = 0;

  CPPUNIT_ASSERT(tract._parameters);
  PetscSection valueSection = tract._parameters->get("value").petscSection();
  Vec          valueVec     = tract._parameters->get("value").localVector();
  PetscScalar *valueArray;
  CPPUNIT_ASSERT(valueSection);CPPUNIT_ASSERT(valueVec);

  err = VecGetArray(valueVec, &valueArray);PYLITH_CHECK_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iPoint) {
    PetscInt vdof, voff;

    err = PetscSectionGetDof(valueSection, v, &vdof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(valueSection, v, &voff);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, vdof);
    
    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar valueE = tractionE[iPoint*spaceDim+d];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, valueArray[voff+d], tolerance);
    } // for
  } // for
  err = VecRestoreArray(valueVec, &valueArray);PYLITH_CHECK_ERROR(err);
} // testCalculate

// ----------------------------------------------------------------------
// Test parameterFields() using 2-D mesh.
void
pylith::faults::TestTractPerturbation::testParameterFields(void)
{ // testParameterFields
  const int spaceDim  = 2;
  const int fiberDimE = 7;
  const PylithScalar parametersE[2*fiberDimE] = {
    0.0, 0.0,   2.0, -1.0,   -1.0, 0.5, 1.5,
    0.0, 0.0,   2.1, -1.1,   -0.8, 0.7, 2.5,
  };

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  TractPerturbation tract;
  _initialize(&mesh, &faultMesh, &tract);
  
  const topology::Fields<topology::Field<topology::SubMesh> >* parameters = tract.parameterFields();

  DM             faultDMMesh = faultMesh.dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  CPPUNIT_ASSERT(faultDMMesh);
  err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);

  CPPUNIT_ASSERT(tract._parameters);
  PetscSection valueSection = tract._parameters->get("value").petscSection();
  Vec          valueVec     = tract._parameters->get("value").localVector();
  PetscScalar *valueArray;
  CPPUNIT_ASSERT(valueSection);CPPUNIT_ASSERT(valueVec);
  PetscSection initialSection = tract._parameters->get("initial").petscSection();
  Vec          initialVec     = tract._parameters->get("initial").localVector();
  PetscScalar *initialArray;
  CPPUNIT_ASSERT(initialSection);CPPUNIT_ASSERT(initialVec);
  PetscSection changeSection = tract._parameters->get("change").petscSection();
  Vec          changeVec     = tract._parameters->get("change").localVector();
  PetscScalar *changeArray;
  CPPUNIT_ASSERT(changeSection);CPPUNIT_ASSERT(changeVec);
  PetscSection changeTimeSection = tract._parameters->get("change time").petscSection();
  Vec          changeTimeVec     = tract._parameters->get("change time").localVector();
  PetscScalar *changeTimeArray;
  CPPUNIT_ASSERT(changeTimeSection);CPPUNIT_ASSERT(changeTimeVec);

  const PylithScalar tolerance = 1.0e-06;
  int iPoint = 0;
  err = VecGetArray(valueVec, &valueArray);PYLITH_CHECK_ERROR(err);
  err = VecGetArray(initialVec, &initialArray);PYLITH_CHECK_ERROR(err);
  err = VecGetArray(changeVec, &changeArray);PYLITH_CHECK_ERROR(err);
  err = VecGetArray(changeTimeVec, &changeTimeArray);PYLITH_CHECK_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iPoint) {
    PetscInt vdof, voff, e = 0;

    err = PetscSectionGetDof(valueSection, v, &vdof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(valueSection, v, &voff);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, vdof);
    PetscInt idof, ioff;

    err = PetscSectionGetDof(initialSection, v, &idof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(initialSection, v, &ioff);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, idof);
    PetscInt cdof, coff;

    err = PetscSectionGetDof(changeSection, v, &cdof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(changeSection, v, &coff);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, cdof);
    PetscInt ctdof, ctoff;

    err = PetscSectionGetDof(changeTimeSection, v, &ctdof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(changeTimeSection, v, &ctoff);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(1, ctdof);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, vdof+idof+cdof+ctdof);

    for(PetscInt d = 0; d < vdof; ++d, ++e) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[iPoint*fiberDimE+e], valueArray[voff+d], tolerance);
    } // for
    for(PetscInt d = 0; d < idof; ++d, ++e) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[iPoint*fiberDimE+e], initialArray[ioff+d], tolerance);
    } // for
    for(PetscInt d = 0; d < cdof; ++d, ++e) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[iPoint*fiberDimE+e], changeArray[coff+d], tolerance);
    } // for
    for(PetscInt d = 0; d < ctdof; ++d, ++e) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[iPoint*fiberDimE+e], changeTimeArray[ctoff+d], tolerance);
    } // for
  } // for
  err = VecRestoreArray(valueVec, &valueArray);PYLITH_CHECK_ERROR(err);
  err = VecRestoreArray(initialVec, &initialArray);PYLITH_CHECK_ERROR(err);
  err = VecRestoreArray(changeVec, &changeArray);PYLITH_CHECK_ERROR(err);
  err = VecRestoreArray(changeTimeVec, &changeTimeArray);PYLITH_CHECK_ERROR(err);
} // testParameterFields

// ----------------------------------------------------------------------
// Test vertexField() using 2-D mesh.
void
pylith::faults::TestTractPerturbation::testVertexField(void)
{ // testVertexField
  const int fiberDimE = 1;
  const PylithScalar fieldE[2*fiberDimE] = {
    1.5,
    2.5,
  };
  const char* label = "traction_change_start_time";

  topology::Mesh mesh;
  topology::SubMesh faultMesh;
  TractPerturbation tract;
  _initialize(&mesh, &faultMesh, &tract);

  const topology::Field<topology::SubMesh>& field = tract.vertexField(label);
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  PetscScalar *array;
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);

  DM             faultDMMesh = faultMesh.dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  CPPUNIT_ASSERT(faultDMMesh);
  err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);

  const PylithScalar tolerance = 1.0e-06;
  int iPoint = 0;
  err = VecGetArray(vec, &array);PYLITH_CHECK_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iPoint) {
    PetscInt dof, off;

    err = PetscSectionGetDof(section, v, &dof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(section, v, &off);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, dof);

    for(PetscInt d = 0; d < dof; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(fieldE[iPoint*fiberDimE+d], array[off+d], tolerance);
    } // for
  } // for
  err = VecRestoreArray(vec, &array);PYLITH_CHECK_ERROR(err);
} // testVertexField

// ----------------------------------------------------------------------
// Initialize TractPerturbation.
void
pylith::faults::TestTractPerturbation::_initialize(topology::Mesh* mesh,
						   topology::SubMesh* faultMesh,
						   TractPerturbation* tract)
{ // _initialize
  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(faultMesh);
  CPPUNIT_ASSERT(tract);
  PetscErrorCode err;

  const char* meshFilename = "data/tri3.mesh";
  const char* faultLabel = "fault";
  const int faultId = 2;
  const char* initialFilename = "data/tri3_initialtract.spatialdb";
  const char* changeFilename = "data/tri3_changetract.spatialdb";
  const PylithScalar orientationVertex[4] = {
      0.0, -1.0,  -1.0, 0.0,
  };

  meshio::MeshIOAscii meshIO;
  meshIO.filename(meshFilename);
  meshIO.debug(false);
  meshIO.interpolate(false);
  meshIO.read(mesh);

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  const int spaceDim = cs.spaceDim();

  DM       dmMesh = mesh->dmMesh(), faultBoundaryDM;
  PetscInt labelSize;
  err = DMPlexGetStratumSize(dmMesh, faultLabel, 1, &labelSize);PYLITH_CHECK_ERROR(err);

  // Create fault mesh
  PetscInt firstFaultVertex    = 0;
  PetscInt firstLagrangeVertex = labelSize;
  PetscInt firstFaultCell      = labelSize;
  DMLabel  groupField;
  const bool useLagrangeConstraints = true;

  err = DMPlexGetLabel(dmMesh, faultLabel, &groupField);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT(groupField);
  if (useLagrangeConstraints) {
    firstFaultCell += labelSize;
  } // if
  ALE::Obj<SieveFlexMesh> faultBoundary = 0;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CohesiveTopology::createFault(faultMesh, faultBoundary, faultBoundaryDM,
                                *mesh, groupField);
  CohesiveTopology::create(mesh, *faultMesh, faultBoundary, faultBoundaryDM,
                           groupField,
                           faultId,
                           firstFaultVertex, firstLagrangeVertex, firstFaultCell,
                           useLagrangeConstraints);
  err = DMDestroy(&faultBoundaryDM);PYLITH_CHECK_ERROR(err);

  // Need to copy coordinates from mesh to fault mesh since we are
  // using create() instead of createParallel().
  _setupFaultCoordinates(mesh, faultMesh);

  // Setup databases
  spatialdata::spatialdb::SimpleDB dbInitial("initial traction");
  spatialdata::spatialdb::SimpleIOAscii ioInitial;
  ioInitial.filename(initialFilename);
  dbInitial.ioHandler(&ioInitial);
  
  spatialdata::spatialdb::SimpleDB dbChange("traction change");
  spatialdata::spatialdb::SimpleIOAscii ioChange;
  ioChange.filename(changeFilename);
  dbChange.ioHandler(&ioChange);

  // Setup fault orientation
  topology::Field<topology::SubMesh> faultOrientation(*faultMesh);
  faultOrientation.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim*spaceDim);
  faultOrientation.allocate();
  PetscSection orientationSection = faultOrientation.petscSection();
  Vec          orientationVec     = faultOrientation.localVector();
  PetscScalar *orientationArray;
  PetscInt           vStart, vEnd;
  err = DMPlexGetDepthStratum(faultMesh->dmMesh(), 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT(orientationSection);CPPUNIT_ASSERT(orientationVec);
  err = VecGetArray(orientationVec, &orientationArray);PYLITH_CHECK_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt odof, ooff;

    err = PetscSectionGetDof(orientationSection, v, &odof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(orientationSection, v, &ooff);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim*spaceDim, odof);
    for(PetscInt d = 0; d < odof; ++d) {
      orientationArray[ooff+d] = orientationVertex[d];
    }
  } // for
  err = VecRestoreArray(orientationVec, &orientationArray);PYLITH_CHECK_ERROR(err);
  
  spatialdata::units::Nondimensional normalizer;

  // setup TractPerturbation
  tract->dbInitial(&dbInitial);
  tract->dbChange(&dbChange);
  
  tract->initialize(*faultMesh, faultOrientation, normalizer);
} // _initialize

// ----------------------------------------------------------------------
// Setup fault coordinates
void
pylith::faults::TestTractPerturbation::_setupFaultCoordinates(topology::Mesh *mesh, topology::SubMesh *faultMesh)
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
  err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetHeightStratum(faultDMMesh, 1, &ffStart, &ffEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexCreateSubpointIS(faultDMMesh, &subpointIS);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetCoordinateSection(faultDMMesh, &fcoordSection);PYLITH_CHECK_ERROR(err);
  err = PetscSectionSetChart(fcoordSection, vStart, vEnd);PYLITH_CHECK_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    err = PetscSectionSetDof(fcoordSection, v, spaceDim);PYLITH_CHECK_ERROR(err);
  }
  err = PetscSectionSetUp(fcoordSection);PYLITH_CHECK_ERROR(err);
  Vec          coordVec, fcoordVec;
  PetscScalar *coords,  *fcoords;
  PetscInt     coordSize;

  err = PetscSectionGetStorageSize(fcoordSection, &coordSize);PYLITH_CHECK_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);PYLITH_CHECK_ERROR(err);
  err = VecCreate(mesh->comm(), &fcoordVec);PYLITH_CHECK_ERROR(err);
  err = VecSetSizes(fcoordVec, coordSize, PETSC_DETERMINE);PYLITH_CHECK_ERROR(err);
  err = VecSetFromOptions(fcoordVec);PYLITH_CHECK_ERROR(err);
  err = ISGetIndices(subpointIS, &points);PYLITH_CHECK_ERROR(err);
  err = VecGetArray(coordVec, &coords);PYLITH_CHECK_ERROR(err);
  err = VecGetArray(fcoordVec, &fcoords);PYLITH_CHECK_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt off, foff;

    // Notice that subpointMap[] does not account for cohesive cells
    err = PetscSectionGetOffset(coordSection, points[v]+(ffEnd-ffStart), &off);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(fcoordSection, v, &foff);PYLITH_CHECK_ERROR(err);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      fcoords[foff+d] = coords[off+d];
    }
  }
  err = ISRestoreIndices(subpointIS, &points);PYLITH_CHECK_ERROR(err);
  err = ISDestroy(&subpointIS);PYLITH_CHECK_ERROR(err);
  err = VecRestoreArray(coordVec, &coords);PYLITH_CHECK_ERROR(err);
  err = VecRestoreArray(fcoordVec, &fcoords);PYLITH_CHECK_ERROR(err);
  err = DMSetCoordinatesLocal(faultDMMesh, fcoordVec);PYLITH_CHECK_ERROR(err);
  err = VecDestroy(&fcoordVec);PYLITH_CHECK_ERROR(err);
} // _setupFaultCoordinates


// End of file 
