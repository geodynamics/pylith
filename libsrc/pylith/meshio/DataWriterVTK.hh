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

/**
 * @file libsrc/meshio/DataWriterVTK.hh
 *
 * @brief Object for writing finite-element data to VTK file.
 *
 * The PETSc VTK viewer collects all fields and then writes them all
 * at once. This means we need to cache fields locally in order to
 * allow the output manager to reuse fields for dimensionalizing,
 * etc. Other writers do not suffer from this restriction, so we
 * implement this functionality in DataWriterVTK.
 */

#if !defined(pylith_meshio_datawritervtk_hh)
#define pylith_meshio_datawritervtk_hh

// Include directives ---------------------------------------------------
#include "DataWriter.hh" // ISA DataWriter

#include "pylith/topology/topologyfwd.hh" // HOLDSA Field
#include "pylith/utils/petscfwd.h" // HASA PetscDM

// DataWriterVTK --------------------------------------------------------
/// Object for writing finite-element data to VTK file.
class pylith::meshio::DataWriterVTK : public DataWriter {
    friend class TestDataWriterVTKMesh; // unit testing
    friend class TestDataWriterVTKMaterial; // unit testing
    friend class TestDataWriterVTKSubmesh; // unit testing
    friend class TestDataWriterVTKBCMesh; // unit testing
    friend class TestDataWriterVTKFaultMesh; // unit testing
    friend class TestDataWriterVTKPoints; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    DataWriterVTK(void);

    /// Destructor
    ~DataWriterVTK(void);

    /** Make copy of this object.
     *
     * @returns Copy of this.
     */
    DataWriter* clone(void) const;

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set filename for VTK file.
     *
     * @param[in] filename Name of VTK file.
     */
    void filename(const char* filename);

    /** Set time format for time stamp in name of VTK file.
     *
     * @param[in] format C style time format for filename.
     */
    void timeFormat(const char* format);

    /** Set value used to normalize time stamp in name of VTK file.
     *
     * Time stamp is divided by this value (time in seconds).
     *
     * @param[in] value Value (time in seconds) used to normalize time stamp in
     * filename.
     */
    void timeConstant(const PylithScalar value);

    /** Set precision of floating point values in output.
     *
     * @param[in] value Precision for floating point values.
     */
    void precision(const int value);

    /** Prepare for writing files.
     *
     * @param[in] mesh Finite-element mesh.
     * @param[in] isInfo True if only writing info values.
     */
    void open(const topology::Mesh& mesh,
              const bool isInfo);

    /// Close output files.
    void close(void);

    /** Prepare file for data at a new time step.
     *
     * @param[in] t Time stamp for new data
     * @param[in] mesh Finite-element mesh.
     */
    void openTimeStep(const PylithScalar t,
                      const topology::Mesh& mesh);

    /// Cleanup after writing data for a time step.
    void closeTimeStep(void);

    /** Write field over vertices to file.
     *
     * @param[in] t Time associated with field.
     * @param[in] subfield Subfield with basis order 1.
     */
    void writeVertexField(const PylithScalar t,
                          const pylith::meshio::OutputSubfield& field);

    /** Write field over cells to file.
     *
     * @param[in] t Time associated with field.
     * @param[in] subfield Subfield with basis order 0.
     */
    void writeCellField(const PylithScalar t,
                        const pylith::meshio::OutputSubfield& subfield);

    // PRIVATE METHODS //////////////////////////////////////////////////////
private:

    /** Copy constructor.
     *
     * @param[in] w Object to copy.
     */
    DataWriterVTK(const DataWriterVTK& w);

    /** Generate filename for VTK file.
     *
     * @param[in] t Time in seconds.
     */
    std::string _vtkFilename(const PylithScalar t) const;

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    const DataWriterVTK& operator=(const DataWriterVTK&); ///< Not implemented

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    /// Time value (in seconds) used to normalize time stamp.
    PylithScalar _timeConstant;
    int _precision; ///< Precision of floating point values in output.

    std::string _filename; ///< Name of VTK file.
    std::string _timeFormat; ///< C style time format for time stamp.

    PetscViewer _viewer; ///< Output file
    PetscDM _dm; ///< Handle to PETSc DM for mesh

    bool _isOpenTimeStep; ///< true if called openTimeStep().
    bool _wroteVertexHeader; ///< True if wrote header for vertex data.
    bool _wroteCellHeader; ///< True if wrote header for cell data

}; // DataWriterVTK

#include "DataWriterVTK.icc" // inline methods

#endif // pylith_meshio_datawritervtk_hh

// End of file
