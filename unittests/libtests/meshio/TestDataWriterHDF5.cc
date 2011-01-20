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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterHDF5.hh" // Implementation of class methods

#include <string.h> // USES strcmp()
#include <iostream> // USES std::cerr

#include "pylith/meshio/HDF5.hh" // USES HDF5

#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 8
#define PYLITH_HDF5_USE_API_18
#endif

// ----------------------------------------------------------------------
#if 0
namespace pylith {
  namespace meshio {
    namespace _TestDataWriterHDF5 {
#endif
      herr_t checkObject(hid_t id,
			 const char* name,
			 H5O_info_t* info,
			 void* data);
#if 0
    } // _TestDataWriterHDF5
  } // meshio
} // pylith
#endif
herr_t
checkObject(hid_t id, 
						 const char* name, 
						 const H5O_info_t* info,
						 void* data)
{ // checkObject
  CPPUNIT_ASSERT(info);  
  CPPUNIT_ASSERT(data);


  hid_t* file = (hid_t*) data;
  CPPUNIT_ASSERT(H5Iis_valid(*file));
  herr_t err = 0;

  switch (info->type) {
  case H5O_TYPE_GROUP : {
    std::cout << "Examining group '" << name << "'." << std::endl;
#if 1
    hid_t group = H5Gopen2(*file, name, H5P_DEFAULT);
    CPPUNIT_ASSERT(group >= 0);
    err = H5Gclose(group);
    CPPUNIT_ASSERT(err >= 0);
#endif
    break;
  } // group
  case H5O_TYPE_DATASET : {
    std::cout << "Examining dataset '" << name << "'." << std::endl;
#if 1
    hid_t dataset = H5Dopen2(*file, name, H5P_DEFAULT);
    CPPUNIT_ASSERT(dataset >= 0);
    err = H5Dclose(dataset);
    CPPUNIT_ASSERT(err >= 0);
#endif
    break;
  } // dataset
  default :
    CPPUNIT_ASSERT(false);
  } // switch

  return 0;
} // checkObject

// ----------------------------------------------------------------------
// Check HDF5 file against archived file.
void
pylith::meshio::TestDataWriterHDF5::checkFile(const char* filename)
{ // checkFile

  const std::string filenameE = "data/" + std::string(filename);

  std::cout << "Checking file '" << filename << "'." << std::endl;

  herr_t err = 0;

  hid_t fileE = H5Fopen(filenameE.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  CPPUNIT_ASSERT(fileE >= 0);

  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  CPPUNIT_ASSERT(file >= 0);

#if defined(PYLITH_HDF5_USE_API_18)

  // Traverse recursively file with expected values.
  err = H5Ovisit(fileE, H5_INDEX_NAME, H5_ITER_NATIVE, checkObject, (void*) &file);
  CPPUNIT_ASSERT(err >= 0);

#else
#endif

  err = H5Fclose(fileE);
  CPPUNIT_ASSERT(err >= 0);

  err = H5Fclose(file);
  CPPUNIT_ASSERT(err >= 0);

} // checkFile


// End of file 
