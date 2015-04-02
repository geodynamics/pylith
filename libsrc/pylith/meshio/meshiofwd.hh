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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
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
    class ExodusII;
    
    class OutputManager;
    class DataWriter;
    class DataWriterVTK;
    class DataWriterHDF5;
    class DataWriterHDF5Ext;
    class CellFilter;
    class CellFilterAvg;
    class VertexFilter;
    class VertexFilterVecNorm;
    class OutputSolnSubset;
    class OutputSolnPoints;

    class HDF5;
    class Xdmf;

  } // meshio
} // pylith


#endif // pylith_meshio_meshiofwd_hh


// End of file 
