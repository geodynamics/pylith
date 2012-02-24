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

#include "TestNeumann.hh" // Implementation of class methods

#include "pylith/bc/Neumann.hh" // USES Neumann

#include "data/NeumannData.hh" // USES NeumannData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/FieldsNew.hh" // USES FieldsNew
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory

#include "data/NeumannDataQuad4.hh" // USES NeumannDataQuad4
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

#include <stdexcept> // USES std::runtime_erro

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumann );

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveMesh;
typedef pylith::topology::SubMesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::SubMesh::RealUniformSection SubRealUniformSection;

typedef pylith::topology::Field<pylith::topology::SubMesh>::RestrictVisitor RestrictVisitor;

// ----------------------------------------------------------------------
namespace pylith {
  namespace bc {
    namespace _TestNeumann {
      const double pressureScale = 4.0;
      const double lengthScale = 1.0; // Mesh coordinates have scale=1.0
      const double timeScale = 0.5;
      const int ncells = 2;
      const int numQuadPts = 2;
      const int spaceDim = 2;

      const double initial[ncells*numQuadPts*spaceDim] = {
	0.3,  0.4,    0.7,  0.6,
	1.3,  1.4,    1.7,  1.6,
      };
      const double rate[ncells*numQuadPts*spaceDim] = {
	-0.2,  -0.1,   0.4, 0.3,
	-1.2,  -1.1,   1.4, 1.3,
      };
      const double rateTime[ncells*numQuadPts] = {
	0.5,   0.8,
	0.6,   0.9,
      };
      const double change[ncells*numQuadPts*spaceDim] = {
	1.3,  1.4,    1.7,  1.6,
	2.3,  2.4,    2.7,  2.6,
      };
      const double changeTime[ncells*numQuadPts] = {
	2.0,  2.4,
	2.1,  2.5,
      };

      const double tValue = 2.2;
      const double valuesRate[ncells*numQuadPts*spaceDim] = {
	-0.34,  -0.17,  0.56,   0.42,
	-1.92,  -1.76,  1.82,   1.69,
      };
      const double valuesChange[ncells*numQuadPts*spaceDim] = {
	1.3,  1.4,   0.0,  0.0,
	2.3,  2.4,   0.0,  0.0,
      };
      const double valuesChangeTH[ncells*numQuadPts*spaceDim] = {
	1.3*0.98,  1.4*0.98,    0.0,  0.0,
	2.3*0.99,  2.4*0.99,    0.0,  0.0,
      };

      // Check values in section against expected values.
      static
      void _checkValues(const double* valuesE,
			const int fiberDimE,
			const topology::Field<topology::SubMesh>& field);
    } // _TestNeumann
  } // bc
} // pylith

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestNeumann::setUp(void)
{ // setUp
  _data = 0;
  _quadrature = new feassemble::Quadrature<topology::SubMesh>();
  CPPUNIT_ASSERT(0 != _quadrature);
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestNeumann::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestNeumann::testConstructor(void)
{ // testConstructor
  Neumann bc;
} // testConstructor

// ----------------------------------------------------------------------
// Test _getLabel().
void
pylith::bc::TestNeumann::test_getLabel(void)
{ // test_getLabel
  Neumann bc;
  
  const std::string& label = "traction bc";
  bc.label(label.c_str());
  CPPUNIT_ASSERT_EQUAL(label, std::string(bc._getLabel()));
} // test_getLabel

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestNeumann::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  Neumann bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  CPPUNIT_ASSERT(0 != _data);

  const topology::SubMesh& boundaryMesh = *bc._boundaryMesh;
  const ALE::Obj<SieveSubMesh>& submesh = boundaryMesh.sieveMesh();

  // Check boundary mesh
  CPPUNIT_ASSERT(!submesh.isNull());

  const int cellDim = boundaryMesh.dimension();
  const int numCorners = _data->numCorners;
  const int spaceDim = _data->spaceDim;
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = submesh->heightStratum(1);
  const int numVertices = submesh->depthStratum(0)->size();
  const int numCells = cells->size();
  const int boundaryDepth = submesh->depth()-1; // depth of boundary cells

  CPPUNIT_ASSERT_EQUAL(_data->cellDim, cellDim);
  CPPUNIT_ASSERT_EQUAL(_data->numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(_data->numCells, numCells);

  const ALE::Obj<SieveMesh::sieve_type>& sieve = submesh->getSieve();
  ALE::ISieveVisitor::PointRetriever<SieveSubMesh::sieve_type> pV(sieve->getMaxConeSize());
  int dp = 0;
  for(SieveSubMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    const int numCorners = submesh->getNumCellCorners(*c_iter, boundaryDepth);
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, numCorners);

    sieve->cone(*c_iter, pV);
    const SieveSubMesh::point_type *cone = pV.getPoints();
    for(int p = 0; p < pV.getSize(); ++p, ++dp) {
      CPPUNIT_ASSERT_EQUAL(_data->cells[dp], cone[p]);
    }
    pV.clear();
  } // for

  // Check traction values
  const int numQuadPts = _data->numQuadPts;
  const int fiberDim = numQuadPts * spaceDim;
  double_array tractionsCell(fiberDim);
  int index = 0;
  CPPUNIT_ASSERT(0 != bc._parameters);
  const ALE::Obj<SubRealUniformSection>& parametersSection =
    bc._parameters->section();
  CPPUNIT_ASSERT(!parametersSection.isNull());
  const int parametersFiberDim = bc._parameters->fiberDim();
  const int initialIndex = bc._parameters->sectionIndex("initial");
  const int initialFiberDim = bc._parameters->sectionFiberDim("initial");

  const double tolerance = 1.0e-06;
  for(SieveSubMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    CPPUNIT_ASSERT_EQUAL(parametersFiberDim,
			 parametersSection->getFiberDimension(*c_iter));
    const double* parametersCell = parametersSection->restrictPoint(*c_iter);
    CPPUNIT_ASSERT(parametersCell);
    CPPUNIT_ASSERT(initialIndex + initialFiberDim <= parametersFiberDim);
    const double* tractionsCell = &parametersCell[initialIndex];
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
      for (int iDim =0; iDim < spaceDim; ++iDim) {
	const double tractionsCellData = _data->tractionsCell[index];
	CPPUNIT_ASSERT_DOUBLES_EQUAL(tractionsCellData,
				     tractionsCell[iQuad*spaceDim+iDim],
				     tolerance);
	++index;
      } // for
  } // for

} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestNeumann::testIntegrateResidual(void)
{ // testIntegrateResidual
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  Neumann bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  topology::Field<topology::Mesh>& residual = fields.get("residual");
  const double t = 0.0;
  bc.integrateResidual(residual, t, &fields);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CPPUNIT_ASSERT(!sieveMesh->depthStratum(0).isNull());

  const double* valsE = _data->valsResidual;
  const int totalNumVertices = sieveMesh->depthStratum(0)->size();
  const int sizeE = _data->spaceDim * totalNumVertices;

  const ALE::Obj<RealSection>& residualSection = residual.section();
  CPPUNIT_ASSERT(!residualSection.isNull());

  const double* vals = residualSection->restrictSpace();
  const int size = residualSection->sizeWithBC();
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  //residual.view("RESIDUAL");

  const double tolerance = 1.0e-06;
  // std::cout << "computed residuals: " << std::endl;
  for (int i=0; i < size; ++i)
    // std::cout << "  " << vals[i] << std::endl;
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i], vals[i], tolerance);
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test _queryDatabases().
void
pylith::bc::TestNeumann::test_queryDatabases(void)
{ // test_queryDatabases
  _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(0 != _quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc, true);

  spatialdata::spatialdb::SimpleDB dbInitial("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename("data/quad4_traction_initial.spatialdb");
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbRate("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbRateIO;
  dbRateIO.filename("data/quad4_traction_rate.spatialdb");
  dbRate.ioHandler(&dbRateIO);
  dbRate.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbChange("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/quad4_traction_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::TimeHistory th("_TestNeumann _queryDatabases");
  th.filename("data/quad4_traction.timedb");

  bc.dbInitial(&dbInitial);
  bc.dbRate(&dbRate);
  bc.dbChange(&dbChange);
  bc.dbTimeHistory(&th);

  const double pressureScale = _TestNeumann::pressureScale;
  const double timeScale = _TestNeumann::timeScale;
  bc._queryDatabases();

  const double tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(0 != bc._parameters);

  // bc._parameters->view("PARAMETERS"); // DEBUGGING

  // Check initial values.
  const topology::Field<topology::SubMesh>& initial = 
    bc._parameters->get("initial");
  _TestNeumann::_checkValues(_TestNeumann::initial, 
			     numQuadPts*spaceDim, initial);

  // Check rate values.
  const topology::Field<topology::SubMesh>& rate = bc._parameters->get("rate");
  _TestNeumann::_checkValues(_TestNeumann::rate, 
			     numQuadPts*spaceDim, rate);

  // Check rate start time.
  const topology::Field<topology::SubMesh>& rateTime = 
    bc._parameters->get("rate time");
  _TestNeumann::_checkValues(_TestNeumann::rateTime, 
			     numQuadPts, rateTime);

  // Check change values.
  const topology::Field<topology::SubMesh>& change = 
    bc._parameters->get("change");
  _TestNeumann::_checkValues(_TestNeumann::change, 
			     numQuadPts*spaceDim, change);

  // Check change start time.
  const topology::Field<topology::SubMesh>& changeTime = 
    bc._parameters->get("change time");
  _TestNeumann::_checkValues(_TestNeumann::changeTime, 
			     numQuadPts, changeTime);
  th.close();
} // test_queryDatabases

// ----------------------------------------------------------------------
// Test _paramsLocalToGlobal().
void
pylith::bc::TestNeumann::test_paramsLocalToGlobal(void)
{ // test_paramsLocalToGlobal
  _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(0 != _quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc, true);

  spatialdata::spatialdb::SimpleDB dbInitial("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename("data/quad4_traction_initial.spatialdb");
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbRate("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbRateIO;
  dbRateIO.filename("data/quad4_traction_rate.spatialdb");
  dbRate.ioHandler(&dbRateIO);
  dbRate.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbChange("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/quad4_traction_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  bc.dbInitial(&dbInitial);
  bc.dbRate(&dbRate);
  bc.dbChange(&dbChange);

  const double pressureScale = _TestNeumann::pressureScale;
  const double timeScale = _TestNeumann::timeScale;
  bc._queryDatabases();
  const double upDir[3] = { 0.0, 0.0, 1.0 };
  bc._paramsLocalToGlobal(upDir);

  const double tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(0 != bc._parameters);

  // Orientation for quad4 is +x, -y for shear and normal tractions.
  CPPUNIT_ASSERT_EQUAL(2, spaceDim); 
  const int ncells = _TestNeumann::ncells;
  double_array valuesE(ncells*numQuadPts*spaceDim);
  
  // Check initial values.
  for (int i=0; i < valuesE.size(); i+=spaceDim) {
    valuesE[i+0] = _TestNeumann::initial[i+0]; // x
    valuesE[i+1] = -_TestNeumann::initial[i+1]; // y
  } // for
  const topology::Field<topology::SubMesh>& initial = 
    bc._parameters->get("initial");
  _TestNeumann::_checkValues(&valuesE[0], 
			     numQuadPts*spaceDim, initial);

  // Check rate values.
  for (int i=0; i < valuesE.size(); i+=spaceDim) {
    valuesE[i+0] = _TestNeumann::rate[i+0]; // x
    valuesE[i+1] = -_TestNeumann::rate[i+1]; // y
  } // for
  const topology::Field<topology::SubMesh>& rate = bc._parameters->get("rate");
  _TestNeumann::_checkValues(&valuesE[0], 
			     numQuadPts*spaceDim, rate);

  // Check change values.
  for (int i=0; i < valuesE.size(); i+=spaceDim) {
    valuesE[i+0] = _TestNeumann::change[i+0]; // x
    valuesE[i+1] = -_TestNeumann::change[i+1]; // y
  } // for
  const topology::Field<topology::SubMesh>& change = 
    bc._parameters->get("change");
  _TestNeumann::_checkValues(&valuesE[0], 
			     numQuadPts*spaceDim, change);
} // test_paramsLocalToGlobal

// ----------------------------------------------------------------------
// Test _calculateValue() with initial value.
void
pylith::bc::TestNeumann::test_calculateValueInitial(void)
{ // test_calculateValueInitial
  _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(0 != _quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc, true);

  spatialdata::spatialdb::SimpleDB dbInitial("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename("data/quad4_traction_initial.spatialdb");
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  bc.dbInitial(&dbInitial);

  const double timeScale = _TestNeumann::timeScale;
  bc._queryDatabases();
  bc._calculateValue(_TestNeumann::tValue/timeScale);

  const double tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(0 != bc._parameters);
  
  // Check values.
  const topology::Field<topology::SubMesh>& value = 
    bc._parameters->get("value");
  _TestNeumann::_checkValues(_TestNeumann::initial, 
			     numQuadPts*spaceDim, value);
} // test_calculateValueInitial

// ----------------------------------------------------------------------
// Test _calculateValue() with rate.
void
pylith::bc::TestNeumann::test_calculateValueRate(void)
{ // test_calculateValueRate
  _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(0 != _quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc, true);

  spatialdata::spatialdb::SimpleDB dbRate("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbRateIO;
  dbRateIO.filename("data/quad4_traction_rate.spatialdb");
  dbRate.ioHandler(&dbRateIO);
  dbRate.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  bc.dbRate(&dbRate);

  const double timeScale = _TestNeumann::timeScale;
  bc._queryDatabases();
  bc._calculateValue(_TestNeumann::tValue/timeScale);

  const double tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(0 != bc._parameters);
  
  // Check values.
  const topology::Field<topology::SubMesh>& value = 
    bc._parameters->get("value");
  _TestNeumann::_checkValues(_TestNeumann::valuesRate,
			     numQuadPts*spaceDim, value);
} // test_calculateValueRate

// ----------------------------------------------------------------------
// Test _calculateValue() with temporal change.
void
pylith::bc::TestNeumann::test_calculateValueChange(void)
{ // test_calculateValueChange
  _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(0 != _quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc, true);

  spatialdata::spatialdb::SimpleDB dbChange("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/quad4_traction_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  bc.dbChange(&dbChange);

  const double timeScale = _TestNeumann::timeScale;
  bc._queryDatabases();
  bc._calculateValue(_TestNeumann::tValue/timeScale);

  const double tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(0 != bc._parameters);
  
  // Check values.
  const topology::Field<topology::SubMesh>& value = 
    bc._parameters->get("value");
  _TestNeumann::_checkValues(_TestNeumann::valuesChange,
			     numQuadPts*spaceDim, value);
} // test_calculateValueChange

// ----------------------------------------------------------------------
// Test _calculateValue() with temporal change w/time history.
void
pylith::bc::TestNeumann::test_calculateValueChangeTH(void)
{ // test_calculateValueChangeTH
  _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(0 != _quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc, true);


  spatialdata::spatialdb::SimpleDB dbChange("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/quad4_traction_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::TimeHistory th("_TestNeumann _queryDatabases");
  th.filename("data/quad4_traction.timedb");

  bc.dbChange(&dbChange);
  bc.dbTimeHistory(&th);

  const double timeScale = _TestNeumann::timeScale;
  bc._queryDatabases();
  bc._calculateValue(_TestNeumann::tValue/timeScale);

  const double tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(0 != bc._parameters);
  
  // Check values.
  const topology::Field<topology::SubMesh>& value = 
    bc._parameters->get("value");
  _TestNeumann::_checkValues(_TestNeumann::valuesChangeTH,
			     numQuadPts*spaceDim, value);
} // test_calculateValueChangeTH

// ----------------------------------------------------------------------
// Test _calculateValue() with initial, rate, and temporal change w/time history.
void
pylith::bc::TestNeumann::test_calculateValueAll(void)
{ // test_calculateValueAll
  _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(0 != _quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc, true);


  spatialdata::spatialdb::SimpleDB dbInitial("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename("data/quad4_traction_initial.spatialdb");
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbRate("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbRateIO;
  dbRateIO.filename("data/quad4_traction_rate.spatialdb");
  dbRate.ioHandler(&dbRateIO);
  dbRate.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbChange("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/quad4_traction_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::TimeHistory th("_TestNeumann _queryDatabases");
  th.filename("data/quad4_traction.timedb");

  bc.dbInitial(&dbInitial);
  bc.dbRate(&dbRate);
  bc.dbChange(&dbChange);
  bc.dbTimeHistory(&th);

  const double timeScale = _TestNeumann::timeScale;
  bc._queryDatabases();
  bc._calculateValue(_TestNeumann::tValue/timeScale);

  const double tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(0 != bc._parameters);
  
  // Check values.
  const int ncells = _TestNeumann::ncells;
  double_array valuesE(ncells*numQuadPts*spaceDim);
  for (int i=0; i < valuesE.size(); ++i)
    valuesE[i] = 
      _TestNeumann::initial[i] +
      _TestNeumann::valuesRate[i] +
      _TestNeumann::valuesChangeTH[i];
  
  const topology::Field<topology::SubMesh>& value = 
    bc._parameters->get("value");
  _TestNeumann::_checkValues(&valuesE[0], numQuadPts*spaceDim, value);
} // test_calculateValueAll

// ----------------------------------------------------------------------
void
pylith::bc::TestNeumann::_preinitialize(topology::Mesh* mesh,
					    Neumann* const bc,
					    const bool useScales) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != _quadrature);
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != bc);

  try {
    // Set up mesh
    meshio::MeshIOAscii iohandler;
    iohandler.filename(_data->meshFilename);
    iohandler.read(mesh);

    // Set up coordinates
    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(mesh->dimension());
    cs.initialize();

    spatialdata::units::Nondimensional normalizer;
    if (useScales) {
      normalizer.lengthScale(_TestNeumann::lengthScale);
      normalizer.pressureScale(_TestNeumann::pressureScale);
      normalizer.timeScale(_TestNeumann::timeScale);
    } // if

    mesh->coordsys(&cs);
    mesh->nondimensionalize(normalizer);

    // Set up quadrature
    _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			    _data->basisDerivRef, _data->numQuadPts, 
			    _data->numBasis, _data->cellDim,
			    _data->quadPts, _data->numQuadPts, _data->cellDim,
			    _data->quadWts, _data->numQuadPts,
			    _data->spaceDim);

    bc->quadrature(_quadrature);
    bc->label(_data->label);
    bc->normalizer(normalizer);
    bc->createSubMesh(*mesh);
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _preinitialize

// ----------------------------------------------------------------------
void
pylith::bc::TestNeumann::_initialize(topology::Mesh* mesh,
				     Neumann* const bc,
				     topology::SolutionFields* fields) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != bc);
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _quadrature);

  try {
    _preinitialize(mesh, bc);

    // Set up database
    spatialdata::spatialdb::SimpleDB db("TestNeumann");
    spatialdata::spatialdb::SimpleIOAscii dbIO;
    dbIO.filename(_data->spatialDBFilename);
    db.ioHandler(&dbIO);
    db.queryType(spatialdata::spatialdb::SimpleDB::LINEAR);

    const double upDir[] = { 0.0, 0.0, 1.0 };

    bc->dbInitial(&db);
    bc->initialize(*mesh, upDir);

    // Set up fields
    CPPUNIT_ASSERT(0 != fields);
    fields->add("residual", "residual");
    fields->add("disp(t), bc(t+dt)", "displacement");
    fields->solutionName("disp(t), bc(t+dt)");

    topology::Field<topology::Mesh>& residual = fields->get("residual");
    residual.newSection(topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
    residual.allocate();
    residual.zero();

    fields->copyLayout("residual");
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _initialize

// ----------------------------------------------------------------------
// Check values in section against expected values.
void
pylith::bc::_TestNeumann::_checkValues(const double* valuesE,
					   const int fiberDimE,
					   const topology::Field<topology::SubMesh>& field)
{ // _checkValues
  assert(0 != valuesE);

  const topology::SubMesh& boundaryMesh = field.mesh();
  const ALE::Obj<SieveSubMesh>& submesh = boundaryMesh.sieveMesh();
  CPPUNIT_ASSERT(!submesh.isNull());
  const ALE::Obj<RealSection>& section = field.section();
  CPPUNIT_ASSERT(!section.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    submesh->heightStratum(1);

  const double scale = field.scale();

  const size_t ncells = _TestNeumann::ncells;
  CPPUNIT_ASSERT_EQUAL(ncells, cells->size());

  // Check values associated with BC.
  int icell = 0;
  const double tolerance = 1.0e-06;
  for(SieveSubMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter, ++icell) {

    const int fiberDim = section->getFiberDimension(*c_iter);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
    
    const double* values = section->restrictPoint(*c_iter);
    for (int iDim=0; iDim < fiberDimE; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[icell*fiberDimE+iDim]/scale,
				   values[iDim], tolerance);
  } // for
} // _checkValues


// End of file 
