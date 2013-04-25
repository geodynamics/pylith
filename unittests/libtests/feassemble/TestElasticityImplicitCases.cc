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

#include "TestElasticityImplicitCases.hh" // Implementation of class methods


#include "pylith/topology/Mesh.hh" // USES Quadrature<Mesh>
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

#include "pylith/feassemble/GeometryLine1D.hh" // USES GeometryLine1D
#include "pylith/feassemble/GeometryTri2D.hh" // USES GeometryTri2D
#include "pylith/feassemble/GeometryTet3D.hh" // USES GeometryTet3D

#include "pylith/materials/ElasticStrain1D.hh" // USES ElasticStrain1D
#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain
#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField


// ----------------------------------------------------------------------
#include "data/ElasticityImplicitData1DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicit1DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicit1DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitData1DLinear();
  _gravityField = 0;
  GeometryLine1D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);
  _material = new materials::ElasticStrain1D;

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticStrain1D"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityImplicitData1DQuadratic.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicit1DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicit1DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitData1DQuadratic();
  _gravityField = 0;
  GeometryLine1D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);
  _material = new materials::ElasticStrain1D;

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticStrain1D"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityImplicitData2DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicit2DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicit2DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitData2DLinear();
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
#include "data/ElasticityImplicitData2DQuadratic.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicit2DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicit2DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitData2DQuadratic();
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
#include "data/ElasticityImplicitData3DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicit3DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicit3DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitData3DLinear();
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
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicit3DQuadratic );
#include "data/ElasticityImplicitData3DQuadratic.hh"

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicit3DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitData3DQuadratic();
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
#include "data/ElasticityImplicitGravData1DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitGrav1DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitGrav1DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitGravData1DLinear();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_gravityField);
  GeometryLine1D geometry;
  _quadrature->refGeometry(&geometry);

  const PylithScalar accScale = _data->lengthScale / (_data->timeScale * _data->timeScale);
  const PylithScalar g = 1.0e8 / accScale;
  const PylithScalar gravityDir[3] = { -1.0, 0.0, 0.0};
  _gravityField->gravAcceleration(g);
  _gravityField->gravityDir(gravityDir[0], gravityDir[1], gravityDir[2]);

  _material = new materials::ElasticStrain1D;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticStrain1D"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityImplicitGravData1DQuadratic.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitGrav1DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitGrav1DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitGravData1DQuadratic();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_gravityField);
  GeometryLine1D geometry;
  _quadrature->refGeometry(&geometry);

  const PylithScalar accScale = _data->lengthScale / (_data->timeScale * _data->timeScale);
  const PylithScalar g = 1.0e8 / accScale;
  const PylithScalar gravityDir[3] = { -1.0, 0.0, 0.0};
  _gravityField->gravAcceleration(g);
  _gravityField->gravityDir(gravityDir[0], gravityDir[1], gravityDir[2]);

  _material = new materials::ElasticStrain1D;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticStrain1D"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityImplicitGravData2DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitGrav2DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitGrav2DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitGravData2DLinear();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_gravityField);
  GeometryTri2D geometry;
  _quadrature->refGeometry(&geometry);

  const PylithScalar accScale = _data->lengthScale / (_data->timeScale * _data->timeScale);
  const PylithScalar g = 1.0e8 / accScale;
  const PylithScalar gravityDir[3] = { 0.0, -1.0, 0.0 };
  _gravityField->gravAcceleration(g);
  _gravityField->gravityDir(gravityDir[0], gravityDir[1], gravityDir[2]);

  _material = new materials::ElasticPlaneStrain;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityImplicitGravData2DQuadratic.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitGrav2DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitGrav2DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitGravData2DQuadratic();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_gravityField);
  GeometryTri2D geometry;
  _quadrature->refGeometry(&geometry);

  const PylithScalar accScale = _data->lengthScale / (_data->timeScale * _data->timeScale);
  const PylithScalar g = 1.0e8 / accScale;
  const PylithScalar gravityDir[3] = { 0.0, -1.0, 0.0 };
  _gravityField->gravAcceleration(g);
  _gravityField->gravityDir(gravityDir[0], gravityDir[1], gravityDir[2]);

  _material = new materials::ElasticPlaneStrain;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityImplicitGravData3DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitGrav3DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitGrav3DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitGravData3DLinear();
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
#include "data/ElasticityImplicitGravData3DQuadratic.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitGrav3DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitGrav3DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitGravData3DQuadratic();
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
