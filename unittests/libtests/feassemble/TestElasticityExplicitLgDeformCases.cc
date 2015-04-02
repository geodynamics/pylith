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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestElasticityExplicitLgDeformCases.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Quadrature<Mesh>
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

#include "pylith/feassemble/GeometryTri2D.hh" // USES GeometryTri2D
#include "pylith/feassemble/GeometryTet3D.hh" // USES GeometryTet3D

#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain
#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

// ----------------------------------------------------------------------
#include "data/ElasticityExplicitLgDeformData2DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitLgDeform2DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitLgDeform2DLinear::setUp(void)
{ // setUp
  TestElasticityExplicitLgDeform::setUp();

  _data = new ElasticityExplicitLgDeformData2DLinear();
  _gravityField = 0;
  CPPUNIT_ASSERT(_quadrature);
  GeometryTri2D geometry;
  _quadrature->refGeometry(&geometry);

  _material = new materials::ElasticPlaneStrain;
  CPPUNIT_ASSERT(_material);

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"), std::string(_data->matType));
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitLgDeformData2DQuadratic.hh"
// :TODO: Update after removing FIAT
//CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitLgDeform2DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitLgDeform2DQuadratic::setUp(void)
{ // setUp
  TestElasticityExplicitLgDeform::setUp();

  _data = new ElasticityExplicitLgDeformData2DQuadratic();
  _gravityField = 0;
  CPPUNIT_ASSERT(_quadrature);
  GeometryTri2D geometry;
  _quadrature->refGeometry(&geometry);

  _material = new materials::ElasticPlaneStrain;
  CPPUNIT_ASSERT(_material);

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"), std::string(_data->matType));
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitLgDeformData3DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitLgDeform3DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitLgDeform3DLinear::setUp(void)
{ // setUp
  TestElasticityExplicitLgDeform::setUp();

  _data = new ElasticityExplicitLgDeformData3DLinear();
  _gravityField = 0;
  CPPUNIT_ASSERT(_quadrature);
  GeometryTet3D geometry;
  _quadrature->refGeometry(&geometry);

  _material = new materials::ElasticIsotropic3D;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticIsotropic3D"), std::string(_data->matType));
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitLgDeformData3DQuadratic.hh"
// :TODO: Update after removing FIAT
//CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitLgDeform3DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitLgDeform3DQuadratic::setUp(void)
{ // setUp
  TestElasticityExplicitLgDeform::setUp();

  _data = new ElasticityExplicitLgDeformData3DQuadratic();
  _gravityField = 0;
  CPPUNIT_ASSERT(_quadrature);
  GeometryTet3D geometry;
  _quadrature->refGeometry(&geometry);

  _material = new materials::ElasticIsotropic3D;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticIsotropic3D"), std::string(_data->matType));
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitLgDeformGravData2DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitLgDeformGrav2DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitLgDeformGrav2DLinear::setUp(void)
{ // setUp
  TestElasticityExplicitLgDeform::setUp();

  _data = new ElasticityExplicitLgDeformGravData2DLinear();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_gravityField);
  GeometryTri2D geometry;
  _quadrature->refGeometry(&geometry);

  const PylithScalar accScale = _data->lengthScale / (_data->timeScale * _data->timeScale);
  const PylithScalar g = 1.0e8 / accScale;
  const PylithScalar gravityDir[] = { 0.0, -1.0, 0.0 };
  _gravityField->gravAcceleration(g);
  _gravityField->gravityDir(gravityDir[0], gravityDir[1], gravityDir[2]);

  _material = new materials::ElasticPlaneStrain;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"), std::string(_data->matType));
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitLgDeformGravData2DQuadratic.hh"
// :TODO: Update after removing FIAT
//CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitLgDeformGrav2DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitLgDeformGrav2DQuadratic::setUp(void)
{ // setUp
  TestElasticityExplicitLgDeform::setUp();

  _data = new ElasticityExplicitLgDeformGravData2DQuadratic();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_gravityField);
  GeometryTri2D geometry;
  _quadrature->refGeometry(&geometry);

  const PylithScalar accScale = _data->lengthScale / (_data->timeScale * _data->timeScale);
  const PylithScalar g = 1.0e8 / accScale;
  const PylithScalar gravityDir[] = { 0.0, -1.0, 0.0 };
  _gravityField->gravAcceleration(g);
  _gravityField->gravityDir(gravityDir[0], gravityDir[1], gravityDir[2]);

  _material = new materials::ElasticPlaneStrain;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"), std::string(_data->matType));
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitLgDeformGravData3DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitLgDeformGrav3DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitLgDeformGrav3DLinear::setUp(void)
{ // setUp
  TestElasticityExplicitLgDeform::setUp();

  _data = new ElasticityExplicitLgDeformGravData3DLinear();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_gravityField);
  GeometryTet3D geometry;
  _quadrature->refGeometry(&geometry);

  const PylithScalar accScale = _data->lengthScale / (_data->timeScale * _data->timeScale);
  const PylithScalar g = 1.0e8 / accScale;
  _gravityField->gravAcceleration(g);

  _material = new materials::ElasticIsotropic3D;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticIsotropic3D"), std::string(_data->matType));
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitLgDeformGravData3DQuadratic.hh"
// :TODO: Update after removing FIAT
//CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitLgDeformGrav3DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitLgDeformGrav3DQuadratic::setUp(void)
{ // setUp
  TestElasticityExplicitLgDeform::setUp();

  _data = new ElasticityExplicitLgDeformGravData3DQuadratic();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_gravityField);
  GeometryTet3D geometry;
  _quadrature->refGeometry(&geometry);

  const PylithScalar accScale = _data->lengthScale / (_data->timeScale * _data->timeScale);
  const PylithScalar g = 1.0e8 / accScale;
  _gravityField->gravAcceleration(g);

  _material = new materials::ElasticIsotropic3D;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticIsotropic3D"), std::string(_data->matType));
} // setUp


// End of file 
