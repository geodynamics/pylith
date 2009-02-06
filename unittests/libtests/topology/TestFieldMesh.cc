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

#include "TestFieldMesh.hh" // Implementation of class methods
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/array.hh" // USES double_array

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldMesh );

// ----------------------------------------------------------------------
namespace pylith {
  namespace topology {
    namespace _TestFieldMesh {
      const int cellDim = 2;
      const int nvertices = 4;
      const int ncells = 1;
      const int ncorners = 4;
      const int cells[] = { 0, 1, 2, 3 };
      const double coordinates[] = {
	0.0, 0.0,
	1.0, 0.0,
	0.0, 1.0,
	1.0, 1.0,
      };
    } // _TestFieldMesh
  } // topology
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldMesh::testConstructor(void)
{ // testConstructor
  Mesh mesh;
  Field<Mesh> field(mesh);
} // testConstructor

// ----------------------------------------------------------------------
// Test section().
void
pylith::topology::TestFieldMesh::testSection(void)
{ // testSection
  Mesh mesh;
  Field<Mesh> field(mesh);

  mesh.createSieveMesh();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(section.isNull());
} // testSection

// ----------------------------------------------------------------------
// Test mesh().
void
pylith::topology::TestFieldMesh::testMesh(void)
{ // testMesh
  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);

  const Mesh& mesh2 = field.mesh();
  CPPUNIT_ASSERT_EQUAL(_TestFieldMesh::cellDim, mesh2.dimension());  
} // testMesh

// ----------------------------------------------------------------------
// Test spaceDim().
void
pylith::topology::TestFieldMesh::testSpaceDim(void)
{ // testSpaceDim
  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);

  CPPUNIT_ASSERT_EQUAL(_TestFieldMesh::cellDim, field.spaceDim());
} // testSpaceDim

// ----------------------------------------------------------------------
// Test newSection().
void
pylith::topology::TestFieldMesh::testNewSection(void)
{ // testNewSection
  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();

  Field<Mesh> field(mesh);
  field.newSection();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
} // testNewSection

// ----------------------------------------------------------------------
// Test newSection(points).
void
pylith::topology::TestFieldMesh::testNewSectionPoints(void)
{ // testNewSectionPoints
  const int fiberDim = 2;

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();

  Field<Mesh> field(mesh);
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  field.newSection(vertices, fiberDim);
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());

  CPPUNIT_ASSERT(!vertices.isNull());
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
} // testNewSectionPoints

// ----------------------------------------------------------------------
// Test newSection(domain).
void
pylith::topology::TestFieldMesh::testNewSectionDomain(void)
{ // testNewSectionDomain
  const int fiberDim = 2;

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);

  const ALE::Obj<Mesh::RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
} // testNewSectionDomain

// ----------------------------------------------------------------------
// Test newSection(chart).
void
pylith::topology::TestFieldMesh::testNewSectionChart(void)
{ // testNewSectionChart
  const int fiberDim = 3;
  const int nconstraints[] = { 0, 2, 1, 3 };
  const int constraints[] = {
              // 0
    0, 3,     // 1
    2,        // 2
    0, 1, 2,  // 3
  };
    
  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();

  // Create field with atlas to use to create new field
  Field<Mesh> fieldSrc(mesh);
  fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  const ALE::Obj<Mesh::RealSection>& sectionSrc = fieldSrc.section();
  CPPUNIT_ASSERT(!sectionSrc.isNull());
  const Mesh::RealSection::chart_type& chart = sectionSrc->getChart();

  Field<Mesh> field(mesh);
  field.newSection(chart, fiberDim);
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
} // testNewSectionChart

// ----------------------------------------------------------------------
// Test newSection(field).
void
pylith::topology::TestFieldMesh::testNewSectionField(void)
{ // testNewSectionField
  const int fiberDim = 3;
  const int nconstraints[] = { 0, 2, 1, 3 };
  const int constraints[] = {
              // 0
    0, 3,     // 1
    2,        // 2
    0, 1, 2,  // 3
  };
    
  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();

  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());

  // Create field with atlas to use to create new field
  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    const ALE::Obj<Mesh::RealSection>& section = fieldSrc.section();
    int iV=0;

    CPPUNIT_ASSERT(!vertices.isNull());
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter)
      section->addConstraintDimension(*v_iter, nconstraints[iV++]);
    fieldSrc.allocate();
  } // Setup source field


  Field<Mesh> field(mesh);
  field.newSection(fieldSrc);
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  int iV = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
    CPPUNIT_ASSERT_EQUAL(nconstraints[iV++], 
			 section->getConstraintDimension(*v_iter));
  } // for
} // testNewSectionField

// ----------------------------------------------------------------------
// Test clear().
void
pylith::topology::TestFieldMesh::testClear(void)
{ // testClear
  Mesh mesh(_TestFieldMesh::cellDim);
  Field<Mesh> field(mesh);

  field.scale(2.0);
  field.vectorFieldType(Field<Mesh>::TENSOR);
  field.addDimensionOkay(true);
  
  field.clear();

  CPPUNIT_ASSERT_EQUAL(1.0, field._scale);
  CPPUNIT_ASSERT_EQUAL(Field<Mesh>::OTHER, field._vecFieldType);
  CPPUNIT_ASSERT_EQUAL(false, field._dimensionsOkay);
} // testClear

// ----------------------------------------------------------------------
// Test allocate().
void
pylith::topology::TestFieldMesh::testAllocate(void)
{ // testAllocate
  const int fiberDim = 3;
  const double scale = 2.0;
  const double valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();

  double_array values(fiberDim);
  int i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  const double tolerance = 1.0e-6;
  i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], values[iDim], tolerance);
    } // for
  } // for
} // testAllocate

// ----------------------------------------------------------------------
// Test zero().
void
pylith::topology::TestFieldMesh::testZero(void)
{ // testZero
  const int fiberDim = 3;
  const double scale = 2.0;
  const double valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();

  double_array values(fiberDim);
  int i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  field.zero();

  const double tolerance = 1.0e-6;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, values[iDim], tolerance);
    } // for
  } // for
} // testZero

// ----------------------------------------------------------------------
// Test complete().
void
pylith::topology::TestFieldMesh::testComplete(void)
{ // testComplete
  const int fiberDim = 3;
  const double scale = 2.0;
  const double valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();

  double_array values(fiberDim);
  int i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  field.complete();

  // Expect no change for this serial test
  i = 0;
  const double tolerance = 1.0e-6;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], values[iDim], tolerance);
    } // for
  } // for
} // testComplete

// ----------------------------------------------------------------------
// Test copy().
void
pylith::topology::TestFieldMesh::testCopy(void)
{ // testCopy
  const int fiberDim = 3;
  const double scale = 2.0;
  const double valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    const ALE::Obj<Mesh::RealSection>& section = fieldSrc.section();
    
    double_array values(fiberDim);
    int i = 0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter) {
      for (int iDim=0; iDim < fiberDim; ++iDim)
	values[iDim] = valuesNondim[i++];
      section->updatePoint(*v_iter, &values[0]);
    } // for
  } // Setup source field

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();

  field.copy(fieldSrc);

  int i = 0;
  double_array values(fiberDim);
  const double tolerance = 1.0e-6;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], values[iDim], tolerance);
    } // for
  } // for
} // testCopy

// ----------------------------------------------------------------------
// Test operator+=().
void
pylith::topology::TestFieldMesh::testOperatorAdd(void)
{ // testOperateAdd
  const int fiberDim = 3;
  const double scale = 2.0;
  const double valuesA[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };
  const double valuesB[] = {
    10.1, 20.2, 30.3,
    10.2, 20.3, 30.4,
    10.3, 20.4, 30.5,
    10.4, 20.5, 30.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    const ALE::Obj<Mesh::RealSection>& section = fieldSrc.section();
    
    double_array values(fiberDim);
    int i = 0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter) {
      for (int iDim=0; iDim < fiberDim; ++iDim)
	values[iDim] = valuesA[i++];
      section->updatePoint(*v_iter, &values[0]);
    } // for
  } // Setup source field

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  { // Setup destination field

    double_array values(fiberDim);
    int i = 0;
    for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter) {
      for (int iDim=0; iDim < fiberDim; ++iDim)
	values[iDim] = valuesB[i++];
      section->updatePoint(*v_iter, &values[0]);
    } // for
  } // Setup destination field

  field += fieldSrc;

  int i = 0;
  double_array values(fiberDim);
  const double tolerance = 1.0e-6;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const double valueE = valuesA[i] + valuesB[i];
      ++i;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, values[iDim], tolerance);
    } // for
  } // for
} // testOperateAdd

// ----------------------------------------------------------------------
// Test dimensionalize().
void
pylith::topology::TestFieldMesh::testDimensionalize(void)
{ // testDimensionalize
  const int fiberDim = 3;
  const double scale = 2.0;
  const double valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  double_array values(fiberDim);
  int i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  field.scale(scale);
  field.addDimensionOkay(true);
  field.dimensionalize();

  i = 0;
  const double tolerance = 1.0e-6;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      const double valueE = valuesNondim[i++]*scale;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, values[iDim], tolerance);
    } // for
  } // for

} // testDimensionalize

// ----------------------------------------------------------------------
// Test view().
void
pylith::topology::TestFieldMesh::testView(void)
{ // testView
  const int fiberDim = 3;
  const double scale = 2.0;
  const double valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  const ALE::Obj<Mesh::RealSection>& section = field.section();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  double_array values(fiberDim);
  int i = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  field.view("Testing view");
} // testView

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_buildMesh(Mesh* mesh)
{ // _buildMesh
  assert(0 != mesh);

  mesh->createSieveMesh(_TestFieldMesh::cellDim);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh->sieveMesh();

  ALE::Obj<Mesh::SieveMesh::sieve_type> sieve = 
    new Mesh::SieveMesh::sieve_type(sieveMesh->comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  ALE::Obj<ALE::Mesh::sieve_type> s = 
    new ALE::Mesh::sieve_type(sieve->comm(), sieve->debug());
  
  const int cellDim = _TestFieldMesh::cellDim;
  const int ncells = _TestFieldMesh::ncells;
  const int* cells = _TestFieldMesh::cells;
  const int nvertices = _TestFieldMesh::nvertices;
  const int ncorners = _TestFieldMesh::ncorners;
  const int spaceDim = _TestFieldMesh::cellDim;
  const double* coordinates = _TestFieldMesh::coordinates;
  const bool interpolate = false;
  ALE::SieveBuilder<ALE::Mesh>::buildTopology(s, cellDim, ncells, (int*) cells,
					      nvertices, interpolate, 
					      ncorners);
  std::map<Mesh::SieveMesh::point_type,Mesh::SieveMesh::point_type> renumbering;
  ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
  sieveMesh->setSieve(sieve);
  sieveMesh->stratify();
  ALE::SieveBuilder<Mesh::SieveMesh>::buildCoordinates(sieveMesh, spaceDim, 
						       coordinates);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);
} // _buildMesh


// End of file 
