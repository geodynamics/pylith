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

/** @file libsrc/meshio/meshiofwd.hh
 *
 * @brief Forward declarations for PyLith meshio objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_meshio_meshiofwd_hh)
#define pylith_meshio_meshiofwd_hh

namespace pylith {
  namespace meshio {

    class BinaryIO;

    class MeshIO;
    class MeshBuilder;
    class MeshIOAscii;
    class MeshIOCubit;
    class MeshIOLagrit;
    class GMVFile;
    class GMVFileAscii;
    class GMVFileBinary;
    class PsetFile;
    class PsetFileAscii;
    class PsetFileBinary;
    
    template<typename mesh_type> class OutputManager;
    template<typename mesh_type> class DataWriter;
    template<typename mesh_type> class DataWriterVTK;
    template<typename mesh_type> class CellFilter;
    template<typename mesh_type> class CellFilterAvg;
    template<typename mesh_type> class VertexFilter;
    template<typename mesh_type> class VertexFilterVecNorm;
    class OutputSolnSubset;

    class UCDFaultFile;

  } // meshio
} // pylith


#endif // pylith_meshio_meshiofwd_hh


// End of file 
