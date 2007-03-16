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

#include "TestMaterial.hh" // Implementation of class methods

#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D
#include "data/MaterialData.hh" // USES MaterialData

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestMaterial );

// ----------------------------------------------------------------------
// Test db()
void
pylith::materials::TestMaterial::testDB(void)
{ // testDB
  const char* label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label);
  
  ElasticIsotropic3D material;
  material.db(&db);
  
  CPPUNIT_ASSERT(0 != material._db);
  CPPUNIT_ASSERT(0 == strcmp(label, material._db->label()));
} // testDB

// ----------------------------------------------------------------------
// Test id()
void
pylith::materials::TestMaterial::testID(void)
{ // testID
  const int id = 346;
  ElasticIsotropic3D material;
  material.id(id);
  
  CPPUNIT_ASSERT(id == material.id());
} // testID

// ----------------------------------------------------------------------
// Test label()
void
pylith::materials::TestMaterial::testLabel(void)
{ // testLabel
  const char* label = "the_database";
  ElasticIsotropic3D material;
  material.label(label);
  
  CPPUNIT_ASSERT(0 == strcmp(label, material.label().c_str()));
} // testLabel
    
// ----------------------------------------------------------------------
// Test initialize()
void
pylith::materials::TestMaterial::testInitialize(void)
{ // testInitialize
  CPPUNIT_ASSERT(false);
} // testInitialize

// ----------------------------------------------------------------------
// Test DBToParameters
void
pylith::materials::TestMaterial::_testDBToParameters(Material* material,
						     const MaterialData& data) const
{ // _testDBToParameters
  CPPUNIT_ASSERT(false);
} // _testDBToParameters

// ----------------------------------------------------------------------
// Test _dbValues() and _numDBValues()
void
pylith::materials::TestMaterial::_testDBValues(Material* material,
					       const MaterialData& data) const
{ // _testDBValues
  const int numDBValues = data.numDBValues;

  CPPUNIT_ASSERT(numDBValues == material->_numDBValues());
  char** const dbValuesE = data.dbValues;
  const char** dbValues = material->_dbValues();
  for (int i=0; i < numDBValues; ++i)
    CPPUNIT_ASSERT(0 == strcmp(dbValuesE[i], dbValues[i]));
} // _testDBValues

// ----------------------------------------------------------------------
// Test _numParameters() and _parameterNames()
void
pylith::materials::TestMaterial::_testParameters(Material* material,
						 const MaterialData& data) const
{ // _testParameters
  const int numParameters = data.numParameters;

  CPPUNIT_ASSERT(numParameters == material->_numParameters());
  char** const namesE = data.parameterNames;
  const char** names = material->_parameterNames();
  for (int i=0; i < numParameters; ++i)
    CPPUNIT_ASSERT(0 == strcmp(namesE[i], names[i]));
} // _testParameters


// End of file 
