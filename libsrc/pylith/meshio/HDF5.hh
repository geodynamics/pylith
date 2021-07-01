// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#if !defined(pylith_meshio_hdf5_hh)
#define pylith_meshio_hdf5_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/utils/arrayfwd.hh" // USES string_vector

#include <hdf5.h> // USES hid_t

// HDF5 -----------------------------------------------------------------
/** High-level serial interface for HDF5 operations.
 *
 * @warning Do not use this object for parallel I/O.
 */
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

  /** Check if HDF5 file has group.
   *
   * @param name Full name of group.
   * @returns True if group exists, false otherwise.
   */
  bool hasGroup(const char* name);

  /** Check if HDF5 file has dataset.
   *
   * @param name Full name of dataset.
   * @returns True if dataset exists, false otherwise.
   */
  bool hasDataset(const char* name);

  /** Get dimensions of dataset.
   *
   * @param dims Array of dimensions. [output]
   * @param ndims Number of dimensions. [output]
   * @param parent Full path of parent dataset for attribute.
   * @param name Name of attribute.
   */
  void getDatasetDims(hsize_t** dims,
		      int* ndims,
		      const char* parent,
		      const char* name);

  /** Get names of datasets in group.
   *
   * @param names Names of datasets.
   * @param group Name of parent.
   */
  void getGroupDatasets(string_vector* names,
			const char* parent);

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

  /** Set scalar attribute (used with external handle to HDF5 file,
   * such as PetscHDF5Viewer).
   *
   * @param h5 HDF5 file.
   * @param parent Full path of parent dataset for attribute.
   * @param name Name of attribute.
   * @param value Attribute value.
   * @param datatype Datatype of scalar.
   */
  static
  void writeAttribute(hid_t h5,
		      const char* parent,
		      const char* name,
		      const void* value,
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

  /** Set string attribute (used with external handle to HDF5 file,
   * such as PetscHDF5Viewer).
   *
   * @param h5 HDF5 file.
   * @param parent Full path of parent dataset for attribute.
   * @param name Name of attribute.
   * @param value String value
   */
  static
  void writeAttribute(hid_t h5,
		      const char* parent,
		      const char* name,
		      const char* value);

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

  /** Read string attribute.
   *
   * @param parent Full path of parent dataset for attribute.
   * @param name Name of attribute.
   * @returns value String value
   */
  std::string readAttribute(const char* parent,
			    const char* name);

  /** Create dataset.
   *
   * @param parent Full path of parent group for dataset.
   * @param name Name of dataset.
   * @param maxDims Maximum dimensions of data.
   * @param dimsChunk Dimensions of data chunks.
   * @param ndims Number of dimensions of data.
   * @param datatype Type of data.
   */
  void createDataset(const char* parent,
		     const char* name,
		     const hsize_t* maxDims,
		     const hsize_t* dimsChunk,
		     const int ndims,
		     hid_t datatype);
  
  /** Append chunk to dataset.
   *
   * @param parent Full path of parent group for dataset.
   * @param name Name of dataset.
   * @param data Data.
   * @param dims Current total dimensions of data.
   * @param dimsChunk Dimension of data chunk to write.
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
   * @param dims Dimensions of chunk.
   * @param ndims Number of dimensions of chunk.
   * @param islice Index of data slice.
   * @param datatype Type of data.
   */
  void readDatasetChunk(const char* parent,
			const char* name,
			char** const data,
			hsize_t** const dimsChunk,
			int* const ndims,
			const int chunk,
			hid_t datatype);

  /** Create dataset associated with data stored in a raw external
   * binary file.
   *
   * @param parent Full path of parent group for dataset.
   * @param name Name of dataset.
   * @param filename Name of external raw data file.
   * @param maxDims Maximum dimensions of data.
   * @param ndims Number of dimensions of data.
   * @param datatype Type of data.
   */
  void createDatasetRawExternal(const char* parent,
				const char* name,
				const char* filename,
				const hsize_t* maxDims,
				const int ndims,
				hid_t datatype);
  
  /** Update the properties of a dataset associated with data stored
   * in a raw external binary file.
   *
   * @param parent Full path of parent group for dataset.
   * @param name Name of dataset.
   * @param dims Dimensions of data.
   * @param ndims Number of dimensions of data.
   */
  void extendDatasetRawExternal(const char* parent,
				const char* name,
				const hsize_t* dims,
				const int ndims);
  
  /** Write dataset comprised of an array of strings. 
   *
   * Data is written as fixed length strings matching the maximum
   * string length.
   *
   * @param parent Full path of parent group for dataset.
   * @param name Name of dataset.
   * @param sarray Array of null terminated C strings.
   * @param nstrings Size of array.
   */
  void writeDataset(const char* parent,
		    const char* name,
		    const char* const* sarray,
		    const int nstrings);

  /** Write dataset comprised of an array of fixed length strings.
   *
   * @param parent Full path of parent group for dataset.
   * @param name Name of dataset.
   * @param sarray Array of null terminated C strings.
   * @param nstrings Size of array.
   * @param slen Fixed string length.
   */
  void writeDataset(const char* parent,
		    const char* name,
		    const char* sarray,
		    const int nstrings,
		    const int slen);

  /** Write dataset comprised of an array of fixed length strings (used with
   * external handle to HDF5 file, such as PetscHDF5Viewer).
   *
   * @param h5 HDF5 file.
   * @param parent Full path of parent group for dataset.
   * @param name Name of dataset.
   * @param sarray Array of null terminated C strings.
   * @param nstrings Size of array.
   * @param slen Fixed string length.
   */
  static
  void writeDataset(hid_t h5,
		    const char* parent,
		    const char* name,
		    const char* sarray,
		    const int nstrings,
		    const int slen);

  /** Read dataset comprised of an array of strings.
   *
   * @param parent Full path of parent group for dataset.
   * @param name Name of dataset.
   * @returns Array of string.
   */
  pylith::string_vector readDataset(const char* parent,
				    const char* name);

// PRIVATE MEMBERS ------------------------------------------------------
private :

  hid_t _file; ///< HDF5 file

}; // HDF5

#endif // pylith_meshio_hdf5_hh

// End of file 
