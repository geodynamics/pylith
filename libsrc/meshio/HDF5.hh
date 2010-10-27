// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#if !defined(pylith_meshio_hdf5_hh)
#define pylith_meshio_hdf5_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

extern "C" {
#include "hdf5.h" // USES hid_t
}

#include <string> // USES std::string

// HDF5 -----------------------------------------------------------------
/// High-level interface for HDF5 operations.
class pylith::meshio::HDF5
{ // HDF5
  friend class TestHDF5; // Unit testing

// PUBLIC METHODS -------------------------------------------------------
public :

  /// Default constructor.
  HDF5(void);

  /** Constructor with filename and mode.
   *
   * @param filename Name of HDF5 file
   * @param mode Mode for HDF5 file
   */
  HDF5(const char* filename,
       hid_t mode);

  /// Destructor
  ~HDF5(void);

  /** Open HDF5 file.
   *
   * @param filename Name of HDF5 file
   * @param mode Mode for HDF5 file
   */
  void open(const char* filename,
	    hid_t mode);

  /// Close HDF5 file.
  void close(void);

  /** Check if HDF5 file is open.
   *
   * @returns True if HDF5 file is open, false otherwise.
   */
  bool isOpen(void) const;

  /** Create group.
   *
   * Create group and leave group open for further operations.
   *
   * @param name Name of group (with absolute path).
   */
  void createGroup(const char* name);

  /** Set scalar attribute.
   *
   * @param parent Full path of parent dataset for attribute.
   * @param name Name of attribute.
   * @param value Attribute value.
   * @param datatype Datatype of scalar.
   */
  void writeAttribute(const char* parent,
		      const char* name,
		      const void* value,
		      hid_t datatype);

  /** Read scalar attribute.
   *
   * @param parent Full path of parent dataset for attribute.
   * @param name Name of attribute.
   * @param datatype Datatype of scalar.
   * @param value Attribute value.
   */
  void readAttribute(const char* parent,
		     const char* name,
		     void* value,
		     hid_t datatype);

  /** Set string attribute.
   *
   * @param parent Full path of parent dataset for attribute.
   * @param name Name of attribute.
   * @param value String value
   */
  void writeAttribute(const char* parent,
		      const char* name,
		      const char* value);

  /** Read string attribute.
   *
   * @param parent Full path of parent dataset for attribute.
   * @param name Name of attribute.
   * @param value String value
   */
  std::string readAttribute(const char* parent,
			    const char* name);

  /** Create dataset.
   *
   * @param parent Full path of parent group for dataset.
   * @param name Name of dataset.
   * @param dims Dimensions of data.
   * @param dimsChunk Dimensions of data chunks.
   * @param ndims Number of dimensions of data.
   * @param datatype Type of data.
   */
  void createDataset(const char* parent,
		     const char* name,
		     const hsize_t* dims,
		     const hsize_t* dimsChunk,
		     const int ndims,
		     hid_t datatype);
  
  /** Create dataset associated with data stored in a raw external
   * binary file.
   *
   * @param parent Full path of parent group for dataset.
   * @param name Name of dataset.
   * @param filename Name of external raw data file.
   * @param dims Dimensions of data.
   * @param ndims Number of dimensions of data.
   * @param datatype Type of data.
   * @returns Dataset identifier.
   */
  void createDatasetRawExternal(const char* parent,
				const char* name,
				const char* filename,
				const hsize_t* dims,
				const int ndims,
				hid_t datatype);
  
  /** Append chunk to dataset.
   *
   * @param parent Full path of parent group for dataset.
   * @param name Name of dataset.
   * @param data Data.
   * @param dims Dimensions of data.
   * @param dimsChunk Dimensions of data chunks.
   * @param ndims Number of dimensions of data.
   * @param chunk Index of data chunk.
   * @param datatype Type of data.
   */
  void writeDatasetChunk(const char* parent,
			 const char* name,
			 const void* data,
			 const hsize_t* dims,
			 const hsize_t* dimsChunk,
			 const int ndims,
			 const int chunk,
			 hid_t datatype);

  /** Read dataset chunk.
   *
   * Currently this method assumes the chunk size (slice along dim=0).
   *
   * @param parent Full path of parent group for dataset.
   * @param name Name of dataset.
   * @param data Data.
   * @param dims Dimensions of data.
   * @param ndims Number of dimensions of data.
   * @param islice Index of data slice.
   * @param datatype Type of data.
   */
  void readDatasetChunk(const char* parent,
			const char* name,
			char** const data,
			hsize_t** const dims,
			int* const ndims,
			const int chunk,
			hid_t datatype);

// PRIVATE MEMBERS ------------------------------------------------------
private :

  hid_t _file; ///< HDF5 file

}; // HDF5

#endif // pylith_meshio_hdf5_hh

// End of file 
