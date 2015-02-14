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

#include "TestElasticityExplicitCases.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Quadrature<Mesh>
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

#include "pylith/feassemble/GeometryTri2D.hh" // USES GeometryTri2D
#include "pylith/feassemble/GeometryTet3D.hh" // USES GeometryTet3D

#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain
#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

// ----------------------------------------------------------------------
#include "data/ElasticityExplicitData2DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicit2DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicit2DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityExplicit::setUp();

  _data = new ElasticityExplicitData2DLinear();
  _gravityField = 0;
  CPPUNIT_ASSERT(_quadrature);
  GeometryTri2D geometry;
  _quadrature->refGeometry(&geometry);

  _material = new materials::ElasticPlaneStrain;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitData2DQuadratic.hh"
// :TODO: Update after removing FIAT
//CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicit2DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicit2DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityExplicit::setUp();

  _data = new ElasticityExplicitData2DQuadratic();
  _gravityField = 0;
  CPPUNIT_ASSERT(_quadrature);
  GeometryTri2D geometry;
  _quadrature->refGeometry(&geometry);

  _material = new materials::ElasticPlaneStrain;
  CPPUNIT_ASSERT(_material);

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitData3DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicit3DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicit3DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityExplicit::setUp();

  _data = new ElasticityExplicitData3DLinear();
  _gravityField = 0;
  CPPUNIT_ASSERT(_quadrature);
  GeometryTet3D geometry;
  _quadrature->refGeometry(&geometry);

  _material = new materials::ElasticIsotropic3D;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticIsotropic3D"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitData3DQuadratic.hh"
// :TODO: Update after removing FIAT
//CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicit3DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicit3DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityExplicit::setUp();

  _data = new ElasticityExplicitData3DQuadratic();
  _gravityField = 0;
  CPPUNIT_ASSERT(_quadrature);
  GeometryTet3D geometry;
  _quadrature->refGeometry(&geometry);

  _material = new materials::ElasticIsotropic3D;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticIsotropic3D"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitGravData2DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitGrav2DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitGrav2DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityExplicit::setUp();

  _data = new ElasticityExplicitGravData2DLinear();
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

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitGravData2DQuadratic.hh"
// :TODO: Update after removing FIAT
//CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitGrav2DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitGrav2DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityExplicit::setUp();

  _data = new ElasticityExplicitGravData2DQuadratic();
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

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitGravData3DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitGrav3DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitGrav3DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityExplicit::setUp();

  _data = new ElasticityExplicitGravData3DLinear();
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

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityExplicitGravData3DQuadratic.hh"
// :TODO: Update after removing FIAT
//CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitGrav3DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitGrav3DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityExplicit::setUp();

  _data = new ElasticityExplicitGravData3DQuadratic();
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

  PYLITH_METHOD_END;
} // setUp


// End of file 
