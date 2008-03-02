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

#if !defined(pylith_meshio_datawritervtkdatamesh_hh)
#define pylith_meshio_datawritervtkdatamesh_hh

namespace pylith {
  namespace meshio {
     class DataWriterVTKDataMesh;
  } // meshio
} // pylith

class pylith::meshio::DataWriterVTKDataMesh
{ // DataWriterVTKDataMesh

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  DataWriterVTKDataMesh(void);

  /// Destructor
  virtual
  ~DataWriterVTKDataMesh(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  char* meshFilename; ///< Name of mesh file.
  char* faultLabel; ///< Name of group of vertices for fault.
  int faultId; ///< Material identifier for fault.
  char* bcLabel; ///< Name of group of vertices for bc.

}; // DataWriterVTKDataMesh

#endif // pylith_meshio_datawritervtkdatamesh_hh


// End of file
