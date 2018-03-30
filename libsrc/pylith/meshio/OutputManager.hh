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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/OutputManager.hh
 *
 * @brief Manager for output of finite-element data.
 */

#if !defined(pylith_meshio_OutputManager_hh)
#define pylith_meshio_OutputManager_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/array.hh" // HASA string_vector

// OutputManager --------------------------------------------------------
/// Manager for output of finite-element data.
class pylith::meshio::OutputManager : public pylith::utils::PyreComponent { // OutputManager
    friend class TestOutputManager;   // unit testing

    // PUBLIC ENUMS /////////////////////////////////////////////////////////
public:

    enum TriggerEnum {
        SKIP_TIMESTEPS=0, ///< Skip X time steps between writes.
        ELAPSED_TIME=1, ///< Skip x time between writes.
    }; // TriggerEnum

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    OutputManager(void);

    /// Destructor
    virtual ~OutputManager(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set trigger for how often to write output.
     *
     * @param[in] flag Flag indicating which method to use for determining how often to write output.
     */
    void trigger(TriggerEnum flag);

    /** Get trigger for how often to write otuput.
     *
     * @returns Flag indicating which method to use for determining how often to write output.
     */
    TriggerEnum trigger(void) const;

    /** Set number of time steps to skip between writes.
     *
     * @param[in] Number of time steps to skip between writes.
     */
    void numTimeStepsSkip(const int value);

    /** Get number of time steps to skip between writes.
     *
     * @returns Number of time steps to skip between writes.
     */
    int numTimeStepsSkip(void) const;

    /** Set elapsed time between writes.
     *
     * @param[in] Elapsed time between writes.
     */
    void timeSkip(const double value);

    /** Get elapsed time between writes.
     *
     * @returns Elapsed time between writes.
     */
    double timeSkip(void) const;

    /** Set writer to write data to file.
     *
     * @param[in] datawriter Writer for data.
     */
    void writer(DataWriter* const datawriter);

    /** Set filter for vertex data.
     *
     * @param[in] filter Filter to apply to vertex data before writing.
     */
    void vertexFilter(VertexFilter* const filter);

    /** Set filter for cell data.
     *
     * @param[in] filter Filter to apply to cell data before writing.
     */
    void cellFilter(CellFilter* const filter);

    /** Set names of vertex information related fields to output.
     *
     * @param[in] names Array of field names.
     * @param[in] numNames Length of array.
     */
    void vertexInfoFields(const char* names[],
                          const int numNames);

    /** Get names of vertex information fields requested for output.
     *
     * @returns Array of field names.
     */
    const pylith::string_vector& vertexInfoFields(void) const;

    /** Set names of vertex data related fields to output.
     *
     * @param[in] names Array of field names.
     * @param[in] numNames Length of array.
     */
    void vertexDataFields(const char* names[],
                          const int numNames);

    /** Get names of vertex data fields requested for output.
     *
     * @returns Array of field names.
     */
    const pylith::string_vector& vertexDataFields(void) const;

    /** Set names of cell information related fields to output.
     *
     * @param[in] names Array of field names.
     * @param[in] numNames Length of array.
     */
    void cellInfoFields(const char* names[],
                        const int numNames);

    /** Get names of cell information fields requested for output.
     *
     * @returns Array of field names.
     */
    const pylith::string_vector& cellInfoFields(void) const;

    /** Set names of cell data related fields to output.
     *
     * @param[in] names Array of field names.
     * @param[in] numNames Length of array.
     */
    void cellDataFields(const char* names[],
                        const int numNames);

    /** Get names of cell data fields requested for output.
     *
     * @returns Array of field names.
     */
    const pylith::string_vector& cellDataFields(void) const;

    /** Verify configuration.
     *
     * @param[in] solution Solution field.
     * @param[in] auxField Auxiliary field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution,
                             const pylith::topology::Field& auxField) const;

    /** Write information.
     *
     * @param[in] auxField Auxiliary field.
     */
    virtual
    void writeInfo(const pylith::topology::Field& auxField);

    /** Write solution at time step.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     * @param[in] auxField Auxiliary field.
     */
    virtual
    void writeTimeStep(const PylithReal t,
                       const PylithInt tindex,
                       const pylith::topology::Field& solution,
                       const pylith::topology::Field& auxField);

    /** Prepare for output.
     *
     * @param[in] mesh Finite-element mesh object.
     * @param[in] isInfo True if only writing info values.
     * @param[in] label Name of label defining cells to include in output (=0 means use all cells in mesh).
     * @param[in] labelId Value of label defining which cells to include.
     */
    virtual
    void open(const topology::Mesh& mesh,
              const bool isInfo,
              const char* label=NULL,
              const int labelId=0);

    /// Close output files.
    virtual
    void close(void);

    /** Setup file for writing fields at time step.
     *
     * @param[in] t Time of time step.
     * @param[in] mesh Finite-element mesh object.
     * @param[in] label Name of label defining cells to include in output (=0 means use all cells in mesh).
     * @param[in] labelId Value of label defining which cells to include.
     */
    virtual
    void openTimeStep(const PylithReal t,
                      const topology::Mesh& mesh,
                      const char* label=NULL,
                      const int labelId=0);

    /// End writing fields at time step.
    virtual
    void closeTimeStep(void);

    /** Append finite-element vertex field to file.
     *
     * @param[in] t Time associated with field.
     * @param[in] field Vertex field.
     * @param[in] mesh Mesh for output.
     */
    virtual
    void appendVertexField(const PylithReal t,
                           topology::Field& field,
                           const topology::Mesh& mesh);

    /** Append finite-element cell field to file.
     *
     * @param[in] t Time associated with field.
     * @param[in] field Cell field.
     * @param[in] label Name of label defining cells to include in output (=0 means use all cells in mesh).
     * @param[in] labelId Value of label defining which cells to include.
     */
    virtual
    void appendCellField(const PylithReal t,
                         topology::Field& field,
                         const char* label=NULL,
                         const int labelId=0);

    /** Check whether we want to write output at time t.
     *
     * @param[in] t Time of proposed write.
     * @param[in] tindex Inxex of current time step.
     * @returns True if output should be written at time t, false otherwise.
     */
    bool shouldWrite(const PylithReal t,
                     const PylithInt tindex);

    /** Get buffer for field.
     *
     * Find the most appropriate buffer that matches field, reusing and reallocating as necessary.
     *
     * @param[in] fieldIn Input field.
     * @param[in] name Name of subfield (optional).
     * @returns Field to use as buffer for outputting field.
     */
    pylith::topology::Field& getBuffer(const pylith::topology::Field& fieldIn,
                                       const char* name=NULL);

    // PROTECTED METHODS ////////////////////////////////////////////////////
protected:

    /** Dimension field.
     *
     * @param[in] fieldIn Field to dimensionalize.
     */
    topology::Field& _dimension(topology::Field& fieldIn);

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    pylith::topology::Fields* _fields;   ///< Container with fields used for output.
    DataWriter* _writer;   ///< Writer for data.
    VertexFilter* _vertexFilter;   ///< Filter applied to vertex data.
    CellFilter* _cellFilter;   ///< Filter applied to cell data.

    pylith::string_vector _vertexInfoFields;
    pylith::string_vector _vertexDataFields;
    pylith::string_vector _cellInfoFields;
    pylith::string_vector _cellDataFields;

    PylithReal _timeSkip; ///< Elapsed time between writes.
    PylithReal _timeWrote; ///< Time when data was previously writtern.
    PylithInt _numTimeStepsSkip; ///< Number of time steps to skip between writes.
    PylithInt _timeStepWrote; ///< Time step when data was previously written.
    TriggerEnum _trigger; ///< Flag indicating whether to use elapsed time or number of time steps when deciding when to write.

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    static const char* _pyreComponent; ///< Name of Pyre component.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    OutputManager(const OutputManager&);   ///< Not implemented.
    const OutputManager& operator=(const OutputManager&);   ///< Not implemented

}; // OutputManager

#endif // pylith_meshio_OutputManager_hh


// End of file
