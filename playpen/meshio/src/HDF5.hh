// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#if !defined(pylith_meshio_hdf5_hh)
#define pylith_meshio_hdf5_hh

namespace pylith {
  namespace meshio {
    class HDF5;
  } // meshio
} // pylith

typedef int hid_t; // HASA hid_t

class pylith::meshio::HDF5
{ // HDF5
  
// PUBLIC METHODS -------------------------------------------------------
public :

  /** Constructor.
   *
   * @param filename Name of HDF5 file
   * @param mode Mode for HDF5 file
   */
  HDF5(const char* filename, hid_t mode);

  /// Destructor
  ~HDF5(void);

  /** Create group.
   *
   * Create group and leave group open for further operations.
   *
   * @param name Name of group (with absolute path).
   * @returns HDF5 group
   */
  hid_t createGroup(const char* name);

  /** Create scalar attribute.
   *
   * @param parent Parent of attribute.
   * @param attrName Name of attribute.
   * @param pValue Pointer to scalar value
   * @param datatype Datatype of scalar.
   */
  void writeAttribute(hid_t parent,
		      const char* name,
		      const void* pValue,
		      hid_t datatype);

  /** Create string attribute.
   *
   * @param parent Parent of attribute.
   * @param attrName Name of attribute.
   * @param value String value
   */
  void writeAttribute(hid_t parent,
		      const char* name,
		      const char* value);

  /** Write dataset.
   *
   * @param parent Parent of dataset.
   * @param name Name of dataset.
   * @param pData Pointer to data.
   * @param dims Dimensions of data.
   * @param ndims Number of dimensions of data.
   * @param datatype Type of data.
   */
  void writeDataset(hid_t parent,
		    const char* name,
		    const void* pData,
		    const int* dims,
		    const int ndims,
		    hid_t datatype);

// PRIVATE MEMBERS ------------------------------------------------------
private :

  hid_t _file; ///< HDF5 file

}; // HDF5

#endif // pylith_meshio_hdf5_hh

// End of file 
