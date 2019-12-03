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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestElasticityImplicitLgDeformCases.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Quadrature<Mesh>
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

#include "pylith/feassemble/GeometryTri2D.hh" // USES GeometryTri2D
#include "pylith/feassemble/GeometryTet3D.hh" // USES GeometryTet3D

#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain
#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField


// ----------------------------------------------------------------------
#include "data/ElasticityImplicitLgDeformData2DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitLgDeform2DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitLgDeform2DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestElasticityImplicitLgDeform::setUp();

  _data = new ElasticityImplicitLgDeformData2DLinear();
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
#include "data/ElasticityImplicitLgDeformData2DQuadratic.hh"
// :TODO: Update after removing FIAT
//CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitLgDeform2DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitLgDeform2DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestElasticityImplicitLgDeform::setUp();

  _data = new ElasticityImplicitLgDeformData2DQuadratic();
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
#include "data/ElasticityImplicitLgDeformData3DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitLgDeform3DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitLgDeform3DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestElasticityImplicitLgDeform::setUp();

  _data = new ElasticityImplicitLgDeformData3DLinear();
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
#include "data/ElasticityImplicitLgDeformData3DQuadratic.hh"
// :TODO: Update after removing FIAT
//CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitLgDeform3DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitLgDeform3DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestElasticityImplicitLgDeform::setUp();

  _data = new ElasticityImplicitLgDeformData3DQuadratic();
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
#include "data/ElasticityImplicitLgDeformGravData2DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitLgDeformGrav2DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitLgDeformGrav2DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestElasticityImplicitLgDeform::setUp();

  _data = new ElasticityImplicitLgDeformGravData2DLinear();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_gravityField);
  GeometryTri2D geometry;
  _quadrature->refGeometry(&geometry);

  const PylithScalar accScale = _data->lengthScale / (_data->timeScale * _data->timeScale);
  const PylithScalar g = 1.0e8 / accScale;
  const PylithScalar gravityDir[3] = { 0.0, -1.0, 0.0 };
  _gravityField->gravityAcc(g);
  _gravityField->gravityDir(gravityDir[0], gravityDir[1], gravityDir[2]);

  _material = new materials::ElasticPlaneStrain;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityImplicitLgDeformGravData2DQuadratic.hh"
// :TODO: Update after removing FIAT
//CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitLgDeformGrav2DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitLgDeformGrav2DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestElasticityImplicitLgDeform::setUp();

  _data = new ElasticityImplicitLgDeformGravData2DQuadratic();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_gravityField);
  GeometryTri2D geometry;
  _quadrature->refGeometry(&geometry);

  const PylithScalar accScale = _data->lengthScale / (_data->timeScale * _data->timeScale);
  const PylithScalar g = 1.0e8 / accScale;
  const PylithScalar gravityDir[3] = { 0.0, -1.0, 0.0 };
  _gravityField->gravityAcc(g);
  _gravityField->gravityDir(gravityDir[0], gravityDir[1], gravityDir[2]);

  _material = new materials::ElasticPlaneStrain;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityImplicitLgDeformGravData3DLinear.hh"
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitLgDeformGrav3DLinear );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitLgDeformGrav3DLinear::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestElasticityImplicitLgDeform::setUp();

  _data = new ElasticityImplicitLgDeformGravData3DLinear();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_gravityField);
  GeometryTet3D geometry;
  _quadrature->refGeometry(&geometry);

  const PylithScalar accScale = _data->lengthScale / (_data->timeScale * _data->timeScale);
  const PylithScalar g = 1.0e8 / accScale;
  _gravityField->gravityAcc(g);

  _material = new materials::ElasticIsotropic3D;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticIsotropic3D"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
#include "data/ElasticityImplicitLgDeformGravData3DQuadratic.hh"
// :TODO: Update after removing FIAT
//CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitLgDeformGrav3DQuadratic );

// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitLgDeformGrav3DQuadratic::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
  
  TestElasticityImplicitLgDeform::setUp();

  _data = new ElasticityImplicitLgDeformGravData3DQuadratic();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_gravityField);
  GeometryTet3D geometry;
  _quadrature->refGeometry(&geometry);

  const PylithScalar accScale = _data->lengthScale / (_data->timeScale * _data->timeScale);
  const PylithScalar g = 1.0e8 / accScale;
  _gravityField->gravityAcc(g);

  _material = new materials::ElasticIsotropic3D;
  CPPUNIT_ASSERT(_material);
  
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticIsotropic3D"), std::string(_data->matType));

  PYLITH_METHOD_END;
} // setUp


// End of file 
