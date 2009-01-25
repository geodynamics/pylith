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

#include "TestField.hh" // Implementation of class methods
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/array.hh" // USES double_array

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestField );

// ----------------------------------------------------------------------
namespace pylith {
  namespace topology {
    namespace _TestField {
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
    } // _TestField
  } // topology
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestField::testConstructor(void)
{ // testConstructor
  Mesh mesh;
  mesh.createSieveMesh();
  Field field(mesh.sieveMesh());

  CPPUNIT_ASSERT(!field._mesh.isNull());
} // testConstructor

// ----------------------------------------------------------------------
// Test section().
void
pylith::topology::TestField::testSection(void)
{ // testSection
  Mesh mesh;
  mesh.createSieveMesh();
  Field field(mesh.sieveMesh());

  const ALE::Obj<SieveMesh::real_section_type>& section = field.section();
  CPPUNIT_ASSERT(section.isNull());

  field.newSection();
  CPPUNIT_ASSERT(!section.isNull());
} // testSection

// ----------------------------------------------------------------------
// Test name().
void 
pylith::topology::TestField::testName(void)
{ // testName
  Mesh mesh;
  mesh.createSieveMesh();
  Field field(mesh.sieveMesh());

  CPPUNIT_ASSERT_EQUAL(std::string("unknown"), std::string(field.name()));

  const std::string name = "field A";
  field.name(name.c_str());
  CPPUNIT_ASSERT_EQUAL(name, std::string(field.name()));
} // testName

// ----------------------------------------------------------------------
// Test vectorFieldType().
void
pylith::topology::TestField::testVectorFieldType(void)
{ // testVectorFieldType
  Mesh mesh;
  mesh.createSieveMesh();
  Field field(mesh.sieveMesh());

  CPPUNIT_ASSERT_EQUAL(Field::OTHER, field.vectorFieldType());

  const Field::VectorFieldEnum ftype = Field::TENSOR;
  field.vectorFieldType(ftype);
  CPPUNIT_ASSERT_EQUAL(ftype, field.vectorFieldType());
} // testVectorFieldType

// ----------------------------------------------------------------------
// Test spaceDim().
void
pylith::topology::TestField::testSpaceDim(void)
{ // testSpaceDim
  Mesh mesh;
  _buildMesh(&mesh);
  Field field(mesh.sieveMesh());

  std::cout << "MESH DIM: " << mesh.sieveMesh()->getDimension() << std::endl;

  CPPUNIT_ASSERT_EQUAL(_TestField::cellDim, field.spaceDim());
} // testSpaceDim

// ----------------------------------------------------------------------
// Test scale().
void
pylith::topology::TestField::testScale(void)
{ // testScale
  Mesh mesh;
  mesh.createSieveMesh();
  Field field(mesh.sieveMesh());

  CPPUNIT_ASSERT_EQUAL(1.0, field.scale());

  const double scale = 4.0;
  field.scale(scale);
  CPPUNIT_ASSERT_EQUAL(scale, field.scale());
} // testScale

// ----------------------------------------------------------------------
// Test addDimensionOkay().
void
pylith::topology::TestField::testAddDimensionOkay(void)
{ // testAddDimensionOkay
  Mesh mesh;
  mesh.createSieveMesh();
  Field field(mesh.sieveMesh());

  CPPUNIT_ASSERT_EQUAL(false, field.addDimensionOkay());

  field.addDimensionOkay(true);
  CPPUNIT_ASSERT_EQUAL(true, field.addDimensionOkay());
} // testAddDimensionOkay

// ----------------------------------------------------------------------
// Test copyLayout().
void
pylith::topology::TestField::testCopyLayout(void)
{ // testCopyLayout
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
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  // Create field with atlas to use to create new field
  Field fieldSrc(sieveMesh);
  { // Setup source field
    fieldSrc.newSection();
    const ALE::Obj<SieveRealSection>& section = fieldSrc.section();
    const int spaceDim = _TestField::cellDim;
    section->setChart(SieveMesh::real_section_type::chart_type(
		  *std::min_element(vertices->begin(), vertices->end()),
		  *std::max_element(vertices->begin(), vertices->end())+1));
    section->setFiberDimension(vertices, fiberDim);
    int iV=0;
    for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter)
      section->addConstraintDimension(*v_iter, nconstraints[iV++]);
    sieveMesh->allocate(section);
  } // Setup source field

  Field field(sieveMesh);
  field.copyLayout(fieldSrc);

  const ALE::Obj<SieveRealSection>& section = field.section();
  int iV = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, section->getFiberDimension(*v_iter));
    CPPUNIT_ASSERT_EQUAL(nconstraints[iV++], 
			 section->getConstraintDimension(*v_iter));
  } // for
} // testCopyLayout

// ----------------------------------------------------------------------
// Test clear().
void
pylith::topology::TestField::testClear(void)
{ // testClear
  Mesh mesh(_TestField::cellDim);
  Field field(mesh.sieveMesh());

  field.scale(2.0);
  field.vectorFieldType(Field::TENSOR);
  field.addDimensionOkay(true);
  
  field.clear();

  CPPUNIT_ASSERT_EQUAL(1.0, field._scale);
  CPPUNIT_ASSERT_EQUAL(Field::OTHER, field._vecFieldType);
  CPPUNIT_ASSERT_EQUAL(false, field._dimensionsOkay);
} // testClear

// ----------------------------------------------------------------------
// Test zero().
void
pylith::topology::TestField::testZero(void)
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
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  Field field(sieveMesh);
  field.newSection();
  const ALE::Obj<SieveRealSection>& section = field.section();
  const int spaceDim = _TestField::cellDim;
    
  section->setChart(SieveMesh::real_section_type::chart_type(
		  *std::min_element(vertices->begin(), vertices->end()),
		  *std::max_element(vertices->begin(), vertices->end())+1));
  section->setFiberDimension(vertices, fiberDim);
  sieveMesh->allocate(section);

  double_array values(fiberDim);
  int i = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  field.zero();

  const double tolerance = 1.0e-6;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
pylith::topology::TestField::testComplete(void)
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
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  Field field(sieveMesh);
  field.newSection();
  const ALE::Obj<SieveRealSection>& section = field.section();
  const int spaceDim = _TestField::cellDim;
    
  section->setChart(SieveMesh::real_section_type::chart_type(
		  *std::min_element(vertices->begin(), vertices->end()),
		  *std::max_element(vertices->begin(), vertices->end())+1));
  section->setFiberDimension(vertices, fiberDim);
  sieveMesh->allocate(section);

  double_array values(fiberDim);
  int i = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
pylith::topology::TestField::testCopy(void)
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
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  Field fieldSrc(sieveMesh);
  { // Setup source field
    fieldSrc.newSection();
    const ALE::Obj<SieveRealSection>& section = fieldSrc.section();
    const int spaceDim = _TestField::cellDim;
    
    section->setChart(SieveMesh::real_section_type::chart_type(
		  *std::min_element(vertices->begin(), vertices->end()),
		  *std::max_element(vertices->begin(), vertices->end())+1));
    section->setFiberDimension(vertices, fiberDim);
    sieveMesh->allocate(section);

    double_array values(fiberDim);
    int i = 0;
    for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter) {
      for (int iDim=0; iDim < fiberDim; ++iDim)
	values[iDim] = valuesNondim[i++];
      section->updatePoint(*v_iter, &values[0]);
    } // for
  } // Setup source field

  Field field(sieveMesh);
  field.newSection();
  const ALE::Obj<SieveRealSection>& section = field.section();
  { // Setup destination field
    const int spaceDim = _TestField::cellDim;
    
    section->setChart(SieveMesh::real_section_type::chart_type(
		  *std::min_element(vertices->begin(), vertices->end()),
		  *std::max_element(vertices->begin(), vertices->end())+1));
    section->setFiberDimension(vertices, fiberDim);
    sieveMesh->allocate(section);
  } // Setup destination field

  field.copy(fieldSrc);

  int i = 0;
  double_array values(fiberDim);
  const double tolerance = 1.0e-6;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
pylith::topology::TestField::testOperatorAdd(void)
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
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  Field fieldSrc(sieveMesh);
  { // Setup source field
    fieldSrc.newSection();
    const ALE::Obj<SieveRealSection>& section = fieldSrc.section();
    const int spaceDim = _TestField::cellDim;
    
    section->setChart(SieveMesh::real_section_type::chart_type(
		  *std::min_element(vertices->begin(), vertices->end()),
		  *std::max_element(vertices->begin(), vertices->end())+1));
    section->setFiberDimension(vertices, fiberDim);
    sieveMesh->allocate(section);

    double_array values(fiberDim);
    int i = 0;
    for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	 v_iter != vertices->end();
	 ++v_iter) {
      for (int iDim=0; iDim < fiberDim; ++iDim)
	values[iDim] = valuesA[i++];
      section->updatePoint(*v_iter, &values[0]);
    } // for
  } // Setup source field

  Field field(sieveMesh);
  field.newSection();
  const ALE::Obj<SieveRealSection>& section = field.section();
  { // Setup destination field
    const int spaceDim = _TestField::cellDim;
    
    section->setChart(SieveMesh::real_section_type::chart_type(
		  *std::min_element(vertices->begin(), vertices->end()),
		  *std::max_element(vertices->begin(), vertices->end())+1));
    section->setFiberDimension(vertices, fiberDim);
    sieveMesh->allocate(section);

    double_array values(fiberDim);
    int i = 0;
    for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
pylith::topology::TestField::testDimensionalize(void)
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
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  Field field(sieveMesh);
  field.newSection();
  const ALE::Obj<SieveRealSection>& section = field.section();
  const int spaceDim = _TestField::cellDim;

  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  section->setChart(SieveMesh::real_section_type::chart_type(
		  *std::min_element(vertices->begin(), vertices->end()),
		  *std::max_element(vertices->begin(), vertices->end())+1));
  section->setFiberDimension(vertices, fiberDim);
  sieveMesh->allocate(section);

  double_array values(fiberDim);
  int i = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
pylith::topology::TestField::testView(void)
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
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);

  Field field(sieveMesh);
  field.newSection();
  const ALE::Obj<SieveRealSection>& section = field.section();
  const int spaceDim = _TestField::cellDim;
    
  section->setChart(SieveMesh::real_section_type::chart_type(
		  *std::min_element(vertices->begin(), vertices->end()),
		  *std::max_element(vertices->begin(), vertices->end())+1));
  section->setFiberDimension(vertices, fiberDim);
  sieveMesh->allocate(section);

  double_array values(fiberDim);
  int i = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
pylith::topology::TestField::_buildMesh(Mesh* mesh)
{ // _buildMesh
  assert(0 != mesh);

  mesh->createSieveMesh(_TestField::cellDim);
  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();

  ALE::Obj<SieveMesh::sieve_type> sieve = 
    new SieveMesh::sieve_type(sieveMesh->comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  ALE::Obj<ALE::Mesh::sieve_type> s = 
    new ALE::Mesh::sieve_type(sieve->comm(), sieve->debug());
  
  const int cellDim = _TestField::cellDim;
  const int ncells = _TestField::ncells;
  const int* cells = _TestField::cells;
  const int nvertices = _TestField::nvertices;
  const int ncorners = _TestField::ncorners;
  const int spaceDim = _TestField::cellDim;
  const double* coordinates = _TestField::coordinates;
  const bool interpolate = false;
  ALE::SieveBuilder<ALE::Mesh>::buildTopology(s, cellDim, ncells, (int*) cells,
					      nvertices, interpolate, 
					      ncorners);
  std::map<SieveMesh::point_type,SieveMesh::point_type> renumbering;
  ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
  sieveMesh->setSieve(sieve);
  sieveMesh->stratify();
  ALE::SieveBuilder<SieveMesh>::buildCoordinates(sieveMesh, spaceDim, 
						 coordinates);

} // _buildMesh

// End of file 
