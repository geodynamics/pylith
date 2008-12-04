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
// Test constructor.
void
pylith::topology::TestField::testConstructor(void)
{ // testConstructor
  Mesh mesh;
  Field field(mesh.sieveMesh());
} // testConstructor

// ----------------------------------------------------------------------
// Test section().
void
pylith::topology::TestField::testSection(void)
{ // testSection
  Mesh mesh;
  Field field(mesh.sieveMesh());

  const ALE::Obj<SieveMesh::real_section_type>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
} // testSection

// ----------------------------------------------------------------------
// Test name().
void 
pylith::topology::TestField::testName(void)
{ // testName
  Mesh mesh;
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
  Field field(mesh.sieveMesh());

  CPPUNIT_ASSERT_EQUAL(0, field.spaceDim());

  const int spaceDim = 2;
  field.spaceDim(spaceDim);
  CPPUNIT_ASSERT_EQUAL(spaceDim, field.spaceDim());
} // testSpaceDim

// ----------------------------------------------------------------------
// Test scale().
void
pylith::topology::TestField::testScale(void)
{ // testScale
  Mesh mesh;
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
  Field field(mesh.sieveMesh());

  CPPUNIT_ASSERT_EQUAL(false, field.addDimensionOkay());

  field.addDimensionOkay(true);
  CPPUNIT_ASSERT_EQUAL(true, field.addDimensionOkay());
} // testAddDimensionOkay

// ----------------------------------------------------------------------
// Test dimensionalize().
void
pylith::topology::TestField::testDimensionalize(void)
{ // testDimensionalize
  Mesh mesh;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  Field field(sieveMesh);

  const int spaceDim = 2;
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

  const int fiberDim = 3;
  const double scale = 2.0;
  const double valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };


  ALE::Obj<SieveMesh::sieve_type> sieve = 
    new SieveMesh::sieve_type(sieveMesh->comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  ALE::Obj<ALE::Mesh::sieve_type> s = 
    new ALE::Mesh::sieve_type(sieve->comm(), sieve->debug());
  
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

  const SieveRealGenSection& section = field.section();
  
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
    section->restrictPoint(*v_iter, &values[0], values.size());
    for (int iDim=0; iDim < fiberDim; ++iDim)
      values[iDim] = valuesNondim[i++];
    section->updatePoint(*v_iter, &values[0]);
  } // for

  field.scale(scale);
  field.spaceDim(spaceDim);
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


// End of file 
