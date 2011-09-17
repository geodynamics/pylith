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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDirichletBC.hh" // Implementation of class methods

#include "pylith/bc/DirichletBC.hh" // USES DirichletBC

#include "data/DirichletData.hh" // USES DirichletData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldsNew.hh" // USES FieldsNew
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/UniformDB.hh" // USES UniformDB

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestDirichletBC );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::RealUniformSection RealUniformSection;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestDirichletBC::setUp(void)
{ // setUp
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestDirichletBC::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestDirichletBC::testConstructor(void)
{ // testConstructor
  DirichletBC bc;
} // testConstructor

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestDirichletBC::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  
  const int numCells = sieveMesh->heightStratum(0)->size();
  const int numFixedDOF = _data->numFixedDOF;
  const size_t numPoints = _data->numConstrainedPts;

  // Check points
  const int offset = numCells;
  if (numFixedDOF > 0) {
    CPPUNIT_ASSERT_EQUAL(numPoints, bc._points.size());
    for (int i=0; i < numPoints; ++i)
      CPPUNIT_ASSERT_EQUAL(_data->constrainedPoints[i]+offset, bc._points[i]);
  } // if

  if (numFixedDOF > 0) {
    // Check values
    CPPUNIT_ASSERT(0 != bc._parameters);
    const ALE::Obj<RealUniformSection>& parametersSection =
      bc._parameters->section();
    CPPUNIT_ASSERT(!parametersSection.isNull());
    const int parametersFiberDim = bc._parameters->fiberDim();
    const int initialIndex = bc._parameters->sectionIndex("initial");
    const int initialFiberDim = bc._parameters->sectionFiberDim("initial");
    CPPUNIT_ASSERT_EQUAL(numFixedDOF, initialFiberDim);

    const PylithScalar tolerance = 1.0e-06;
    for (int i=0; i < numPoints; ++i) {
      const int p_value = _data->constrainedPoints[i]+offset;
      CPPUNIT_ASSERT_EQUAL(parametersFiberDim, 
			   parametersSection->getFiberDimension(p_value));
      const PylithScalar* parametersVertex = 
	parametersSection->restrictPoint(p_value);
      CPPUNIT_ASSERT(parametersVertex);
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) 
	CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->valuesInitial[i*numFixedDOF+iDOF],
				     parametersVertex[initialIndex+iDOF],
				     tolerance);
    } // for
    
    // Check rate of change
    const int rateIndex = bc._parameters->sectionIndex("rate");
    const int rateFiberDim = bc._parameters->sectionFiberDim("rate");
    CPPUNIT_ASSERT_EQUAL(numFixedDOF, rateFiberDim);
    
    for (int i=0; i < numPoints; ++i) {
      const int p_value = _data->constrainedPoints[i]+offset;
      CPPUNIT_ASSERT_EQUAL(parametersFiberDim, 
			   parametersSection->getFiberDimension(p_value));
      const PylithScalar* parametersVertex = 
	parametersSection->restrictPoint(p_value);
      CPPUNIT_ASSERT(parametersVertex);
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) 
	CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->valueRate, 
				     parametersVertex[rateIndex+iDOF],
				     tolerance);
    } // for
  } // if
} // testInitialize

// ----------------------------------------------------------------------
// Test numDimConstrained().
void
pylith::bc::TestDirichletBC::testNumDimConstrained(void)
{ // testNumDimConstrained
  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, bc.numDimConstrained());
} // testNumDimConstrained

// ----------------------------------------------------------------------
// Test setConstraintSizes().
void
pylith::bc::TestDirichletBC::testSetConstraintSizes(void)
{ // testSetConstraintSizes
  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
		 sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  
  const int fiberDim = _data->numDOF;
  const int spaceDim = mesh.dimension();
  topology::Field<topology::Mesh> field(mesh);
  field.newSection(vertices, fiberDim);
  field.splitDefault();
  const ALE::Obj<RealSection>& fieldSection = field.section();
  CPPUNIT_ASSERT(!fieldSection.isNull());

  bc.setConstraintSizes(field);

  const int numCells = sieveMesh->heightStratum(0)->size();
  const int offset = numCells;
  int iConstraint = 0;
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    if (*v_iter != _data->constrainedPoints[iConstraint] + offset) {
      CPPUNIT_ASSERT_EQUAL(_data->numDOF,
			   fieldSection->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(0,
			   fieldSection->getConstraintDimension(*v_iter));
      for (int fibration=0; fibration < spaceDim; ++fibration) {
        CPPUNIT_ASSERT_EQUAL(1,
            fieldSection->getFiberDimension(*v_iter,
                fibration));
        CPPUNIT_ASSERT_EQUAL(0,
            fieldSection->getConstraintDimension(*v_iter,
                fibration));
      } // for
    } else {
      CPPUNIT_ASSERT_EQUAL(_data->numDOF,
			   fieldSection->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, 
			   fieldSection->getConstraintDimension(*v_iter));
      for (int fibration=0; fibration < spaceDim; ++fibration) {
        CPPUNIT_ASSERT_EQUAL(1,
            fieldSection->getFiberDimension(*v_iter,
                fibration));
        bool isConstrained = false;
        for (int iDOF=0; iDOF < _data->numFixedDOF; ++iDOF)
          if (fibration == _data->fixedDOF[iDOF])
            isConstrained = true;
        const int constraintDimE = (!isConstrained) ? 0 : 1;
        CPPUNIT_ASSERT_EQUAL(constraintDimE,
            fieldSection->getConstraintDimension(*v_iter,
                fibration));
      } // for
      ++iConstraint;
    } // if/else
  } // for
} // testSetConstraintSizes

// ----------------------------------------------------------------------
// Test setConstraints().
void
pylith::bc::TestDirichletBC::testSetConstraints(void)
{ // testSetConstraints
  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
		 sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  
  const int spaceDim = mesh.dimension();
  const int fiberDim = _data->numDOF;
  topology::Field<topology::Mesh> field(mesh);
  field.newSection(vertices, fiberDim);
  field.splitDefault();
  const ALE::Obj<RealSection>& fieldSection = field.section();
  CPPUNIT_ASSERT(!fieldSection.isNull());

  bc.setConstraintSizes(field);
  field.allocate();
  bc.setConstraints(field);

  const int numCells = sieveMesh->heightStratum(0)->size();
  const int offset = numCells;
  int iConstraint = 0;
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int* fixedDOF = fieldSection->getConstraintDof(*v_iter);
    if (*v_iter != _data->constrainedPoints[iConstraint] + offset) {
      CPPUNIT_ASSERT_EQUAL(0, fieldSection->getConstraintDimension(*v_iter));
    } else {
      CPPUNIT_ASSERT(0 != fixedDOF);
      CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, 
			   fieldSection->getConstraintDimension(*v_iter));
      for (int iDOF=0; iDOF < _data->numFixedDOF; ++iDOF)
	CPPUNIT_ASSERT_EQUAL(_data->fixedDOF[iDOF], fixedDOF[iDOF]);
      ++iConstraint;
    } // if/else
  } // for

  // Check fibrations for split fields.
  iConstraint = 0;
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    if (*v_iter != _data->constrainedPoints[iConstraint] + offset) {
      for (int fibration=0; fibration < spaceDim; ++fibration)
        CPPUNIT_ASSERT_EQUAL(0,
            fieldSection->getConstraintDimension(*v_iter,
                fibration));
    } else {
      for (int fibration=0; fibration < spaceDim; ++fibration) {
        bool isConstrained = false;
        for (int iDOF=0; iDOF < _data->numFixedDOF; ++iDOF)
          if (fibration == _data->fixedDOF[iDOF])
          isConstrained = true;
        if (isConstrained) {
          CPPUNIT_ASSERT_EQUAL(1,
              fieldSection->getConstraintDimension(*v_iter,
                  fibration));
          const int* fixedDOF = fieldSection->getConstraintDof(*v_iter,
            fibration);
          CPPUNIT_ASSERT(0 != fixedDOF);
          CPPUNIT_ASSERT_EQUAL(0, fixedDOF[0]);
        } else
          CPPUNIT_ASSERT_EQUAL(0,
              fieldSection->getConstraintDimension(*v_iter,
                  fibration));
      } // for
      ++iConstraint;
    } // if/else
  } // for
} // testSetConstraints

// ----------------------------------------------------------------------
// Test setField().
void
pylith::bc::TestDirichletBC::testSetField(void)
{ // testSetField
  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
		 sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  
  const int fiberDim = _data->numDOF;
  topology::Field<topology::Mesh> field(mesh);
  field.newSection(vertices, fiberDim);
  const ALE::Obj<RealSection>& fieldSection = field.section();
  CPPUNIT_ASSERT(!fieldSection.isNull());

  bc.setConstraintSizes(field);
  field.allocate();
  bc.setConstraints(field);

  const PylithScalar tolerance = 1.0e-06;

  // All values should be zero.
  field.zero();
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int fiberDim = fieldSection->getFiberDimension(*v_iter);
    const RealSection::value_type* values = 
      sieveMesh->restrictClosure(fieldSection, *v_iter);
    for (int i=0; i < fiberDim; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[i], tolerance);
  } // for

  // Only unconstrained values should be zero.
  const PylithScalar t = 1.0;
  bc.setField(t, field);

  // Create list of unconstrained DOF at constrained DOF
  const int numFreeDOF = _data->numDOF - _data->numFixedDOF;
  int_array freeDOF(numFreeDOF);
  int index = 0;
  for (int iDOF=0; iDOF < _data->numDOF; ++iDOF) {
    bool free = true;
    for (int iFixed=0; iFixed < _data->numFixedDOF; ++iFixed)
      if (iDOF == _data->fixedDOF[iFixed])
	free = false;
    if (free)
      freeDOF[index] = iDOF;
  } // for

  const int numCells = sieveMesh->heightStratum(0)->size();
  const int offset = numCells;
  const int numFixedDOF = _data->numFixedDOF;
  int iConstraint = 0;
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int fiberDim = fieldSection->getFiberDimension(*v_iter);
    const RealSection::value_type* values = 
      sieveMesh->restrictClosure(fieldSection, *v_iter);

    if (*v_iter != _data->constrainedPoints[iConstraint] + offset) {
      // unconstrained point
      for (int i=0; i < fiberDim; ++i)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[i], tolerance);
    } else {
      // constrained point

      // check unconstrained DOF
      for (int iDOF=0; iDOF < numFreeDOF; ++iDOF)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[freeDOF[iDOF]], tolerance);

      // check constrained DOF
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
	const int index = iConstraint * numFixedDOF + iDOF;
	const PylithScalar valueE = (t > _data->tRef) ?
	  _data->valuesInitial[index] + (t-_data->tRef)*_data->valueRate :
	  _data->valuesInitial[index];
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, values[_data->fixedDOF[iDOF]],
				     tolerance);
      } // for
      ++iConstraint;
    } // if/else
  } // for
} // testSetField

// ----------------------------------------------------------------------
// Test setFieldIncr().
void
pylith::bc::TestDirichletBC::testSetFieldIncr(void)
{ // testSetFieldIncr
  topology::Mesh mesh;
  DirichletBC bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);
  bc.useSolnIncr(true);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
		 sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  
  const int fiberDim = _data->numDOF;
  topology::Field<topology::Mesh> field(mesh);
  field.newSection(vertices, fiberDim);
  const ALE::Obj<RealSection>& fieldSection = field.section();
  CPPUNIT_ASSERT(!fieldSection.isNull());

  bc.setConstraintSizes(field);
  field.allocate();
  bc.setConstraints(field);

  const PylithScalar tolerance = 1.0e-06;

  // All values should be zero.
  field.zero();
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int fiberDim = fieldSection->getFiberDimension(*v_iter);
    const RealSection::value_type* values = 
      sieveMesh->restrictClosure(fieldSection, *v_iter);
    for (int i=0; i < fiberDim; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[i], tolerance);
  } // for

  // Only unconstrained values should be zero.
  const PylithScalar t0 = 1.0;
  const PylithScalar t1 = 2.0;
  bc.setFieldIncr(t0, t1, field);

  // Create list of unconstrained DOF at constrained DOF
  const int numFreeDOF = _data->numDOF - _data->numFixedDOF;
  int_array freeDOF(numFreeDOF);
  int index = 0;
  for (int iDOF=0; iDOF < _data->numDOF; ++iDOF) {
    bool free = true;
    for (int iFixed=0; iFixed < _data->numFixedDOF; ++iFixed)
      if (iDOF == _data->fixedDOF[iFixed])
	free = false;
    if (free)
      freeDOF[index] = iDOF;
  } // for

  const int numCells = sieveMesh->heightStratum(0)->size();
  const int offset = numCells;
  const int numFixedDOF = _data->numFixedDOF;
  int iConstraint = 0;
  for (SieveMesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    const int fiberDim = fieldSection->getFiberDimension(*v_iter);
    const RealSection::value_type* values = 
      sieveMesh->restrictClosure(fieldSection, *v_iter);

    if (*v_iter != _data->constrainedPoints[iConstraint] + offset) {
      // unconstrained point
      for (int i=0; i < fiberDim; ++i)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[i], tolerance);
    } else {
      // constrained point

      // check unconstrained DOF
      for (int iDOF=0; iDOF < numFreeDOF; ++iDOF)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[freeDOF[iDOF]], tolerance);

      // check constrained DOF
      for (int iDOF=0; iDOF < numFixedDOF; ++iDOF) {
	const int index = iConstraint * numFixedDOF + iDOF;
	const PylithScalar valueE = (t0 > _data->tRef) ? (t1-t0)*_data->valueRate :
	  (t1 > _data->tRef) ? (t1-_data->tRef)*_data->valueRate : 0.0;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, values[_data->fixedDOF[iDOF]],
				     tolerance);
      } // for
      ++iConstraint;
    } // if/else
  } // for
} // testSetFieldIncr

// ----------------------------------------------------------------------
void
pylith::bc::TestDirichletBC::_initialize(topology::Mesh* mesh,
					 DirichletBC* const bc) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != bc);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);

  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  mesh->nondimensionalize(normalizer);

  spatialdata::spatialdb::SimpleDB db("TestDirichletBC initial");
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename(_data->dbFilename);
  db.ioHandler(&dbIO);
  db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::UniformDB dbRate("TestDirichletBC rate");
  const int numValues = 4;
  const char* names[] = { 
    "displacement-rate-x", 
    "displacement-rate-y", 
    "displacement-rate-z",
    "rate-start-time"};
  const char* units[] = { 
    "m", 
    "m", 
    "m",
    "s"};
  const double values[numValues] = { 
    _data->valueRate,
    _data->valueRate,
    _data->valueRate,
    _data->tRef,
  };
  dbRate.setData(names, units, values, numValues);

  const double upDir[] = { 0.0, 0.0, 1.0 };

  bc->label(_data->label);
  bc->dbInitial(&db);
  bc->dbRate(&dbRate);
  bc->bcDOF(_data->fixedDOF, _data->numFixedDOF);
  bc->initialize(*mesh, upDir);
} // _initialize


// End of file 
