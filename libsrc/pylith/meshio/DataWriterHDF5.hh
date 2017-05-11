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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/DataWriterHDF5.hh
 *
 * @brief Object for writing finite-element data to HDF5 file.
 *
 * HDF5 schema for PyLith output.
 *
 * / - root group
 *   geometry - group
       coordsys - attribute string with coordinate system
 *     vertices - dataset [nvertices, spacedim]
 *   topology - group
 *     cell_type - attribute string with cell type
 *     cells - dataset [ncells, ncorners]
 *   vertex_fields - group
 *     VERTEX_FIELD (name of vertex field) - dataset
 *       [ntimesteps, nvertices, fiberdim]
 *   cell_fields - group
 *     CELL_FIELD (name of cell field) - dataset
 *       [ntimesteps, ncells, fiberdim]
 *   time - dataset
 *     [ntimesteps]
 *   stations - dataset [optional]
 *     [nvertices, 64]
 */

#if !defined(pylith_meshio_datawriterhdf5_hh)
#define pylith_meshio_datawriterhdf5_hh

// Include directives ---------------------------------------------------
#include "DataWriter.hh" // ISA DataWriter

#include "pylith/utils/petscfwd.h" // HASA PetscVec

#include <string> // USES std::string
#include <map> // HASA std::map

// DataWriterHDF5 --------------------------------------------------------
/// Object for writing finite-element data to HDF5 file.
class pylith::meshio::DataWriterHDF5 : public DataWriter { // DataWriterHDF5
    friend class TestDataWriterHDF5Mesh;   // unit testing
    friend class TestDataWriterHDF5SubMesh;   // unit testing
    friend class TestDataWriterHDF5Points;   // unit testing
    friend class TestDataWriterHDF5BCMesh;   // unit testing
    friend class TestDataWriterHDF5FaultMesh;   // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    DataWriterHDF5(void);

    /// Destructor
    ~DataWriterHDF5(void);

    /** Make copy of this object.
     *
     * @returns Copy of this.
     */
    DataWriter* clone(void) const;

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set filename for HDF5 file.
     *
     * @param[in] filename Name of HDF5 file.
     */
    void filename(const char* filename);

    /** Open output file.
     *
     * @param[in] mesh Finite-element mesh.
     * @param[in] isInfo True if only writing info values.
     * @param[in] label Name of label defining cells to include in output
     *   (=0 means use all cells in mesh).
     * @param[in] labelId Value of label defining which cells to include.
     */
    void open(const topology::Mesh& mesh,
              const bool isInfo,
              const char* label=0,
              const int labelId=0);

    /// Close output files.
    void close(void);

    /** Write field over vertices to file.
     *
     * @param[in] t Time associated with field.
     * @param[in] field Field over vertices.
     * @param[in] mesh Mesh associated with output.
     */
    void writeVertexField(const PylithScalar t,
                          topology::Field& field,
                          const topology::Mesh& mesh);

    /** Write field over cells to file.
     *
     * @param[in] t Time associated with field.
     * @param[in] field Field over cells.
     * @param[in] label Name of label defining cells to include in output
     *   (=0 means use all cells in mesh).
     * @param[in] labelId Value of label defining which cells to include.
     */
    void writeCellField(const PylithScalar t,
                        topology::Field& field,
                        const char* label=0,
                        const int labelId=0);

    /** Write dataset with names of points to file.
     *
     * @param[in] names Array with name for each point, e.g., station name.
     * @param[in] mesh Finite-element mesh.
     *
     * Primarily used with OutputSolnPoints.
     */
    void writePointNames(const pylith::string_vector& names,
                         const topology::Mesh& mesh);

    // PRIVATE METHODS //////////////////////////////////////////////////////
private:

    /** Copy constructor.
     *
     * @param[in] w Object to copy.
     */
    DataWriterHDF5(const DataWriterHDF5& w);

    /// Generate filename for HDF5 file.
    std::string _hdf5Filename(void) const;

    /** Write time stamp to file.
     *
     * @param[in] t Time in seconds.
     * @param[in] commRank Processor rank in MPI communicator.
     */
    void _writeTimeStamp(const PylithScalar t,
                         const int commRank);

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    const DataWriterHDF5& operator=(const DataWriterHDF5&);   ///< Not implemented

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    std::string _filename;   ///< Name of HDF5 file.
    PetscViewer _viewer;   ///< Output file.
    PetscVec _tstamp;   ///< Single value vector holding time stamp.

    std::map<std::string, int> _timesteps;   ///< # of time steps written per field.
    int _tstampIndex;   ///< Index of last time stamp written.

    static const char* _pyreComponent; ///< Name of Pyre component.

}; // DataWriterHDF5

#include "DataWriterHDF5.icc" // inline methods

#endif // pylith_meshio_datawriterhdf5_hh


// End of file
