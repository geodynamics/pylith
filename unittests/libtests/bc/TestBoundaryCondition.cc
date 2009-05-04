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

#include "TestBoundaryCondition.hh" // Implementation of class methods

#include "pylith/bc/DirichletBC.hh" // USES DirichletBC

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB

#include <string> // USES std::string

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryCondition );

// ----------------------------------------------------------------------
// Test label().
void
pylith::bc::TestBoundaryCondition::testLabel(void)
{ // testLabel
  const std::string& label = "the_database";
  DirichletBC bc;
  bc.label(label.c_str());
  
  CPPUNIT_ASSERT_EQUAL(label, std::string(bc.label()));
} // testLabel
    
// ----------------------------------------------------------------------
// Test db().
void
pylith::bc::TestBoundaryCondition::testDB(void)
{ // testDB
  const std::string& label = "my db";
  spatialdata::spatialdb::SimpleDB db(label.c_str());
  DirichletBC bc;
  bc.db(&db);
  
  CPPUNIT_ASSERT(0 != bc._db);
  CPPUNIT_ASSERT_EQUAL(label, std::string(bc._db->label()));
} // testDB
    

// End of file 
