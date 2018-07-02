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

#include "pylith/feassemble/Observer.hh" // ISA Observer
#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/array.hh" // HASA string_vector

// OutputManager --------------------------------------------------------
/// Manager for output of finite-element data.
class pylith::meshio::OutputManager :
    public pylith::utils::PyreComponent,
    public pylith::feassemble::Observer {
    friend class TestOutputManager;   // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    OutputManager(void);

    /// Destructor
    virtual ~OutputManager(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set output trigger for how often to write output.
     *
     * @param[in] otrigger Output trigger.
     */
    void trigger(pylith::meshio::OutputTrigger* const otrigger);

    /** Get trigger for how often to write otuput.
     *
     * @returns Output trigger.
     */
    const pylith::meshio::OutputTrigger* trigger(void) const;

    /** Set writer to write data to file.
     *
     * @param[in] datawriter Writer for data.
     */
    void writer(pylith::meshio::DataWriter* const datawriter);

    /** Set filter for vertex data.
     *
     * @param[in] filter Filter to apply to vertex data before writing.
     */
    void fieldFilter(pylith::meshio::FieldFilter* const filter);

    /** Set names of information fields requested for output.
     *
     * @param[in] names Array of field names.
     * @param[in] numNames Length of array.
     */
    void infoFields(const char* names[],
                    const int numNames);

    /** Get names of information fields requested for output.
     *
     * @returns Array of field names.
     */
    const pylith::string_vector& infoFields(void) const;

    /** Set names of data fields requested for output.
     *
     * @param[in] names Array of field names.
     * @param[in] numNames Length of array.
     */
    void dataFields(const char* names[],
                    const int numNames);

    /** Get names of data fields requested for output.
     *
     * @returns Array of field names.
     */
    const pylith::string_vector& dataFields(void) const;

    /** Receive update from subject.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     * @param[in] infoOnly Flag is true if this update is before solution is available (e.g., after initialization).
     */
    virtual
    void update(const PylithReal t,
                const PylithInt tindex,
                const pylith::topology::Field& solution,
                const bool infoOnly=false);

    // PROTECTED METHODS ////////////////////////////////////////////////////
protected:

    /// Write diagnostic information.
    virtual
    void _writeInfo(void);

    /** Prepare for output.
     *
     * @param[in] mesh Finite-element mesh object.
     * @param[in] isInfo True if only writing info values.
     */
    virtual
    void _open(const pylith::topology::Mesh& mesh,
               const bool isInfo);

    /// Close output files.
    virtual
    void _close(void);

    /** Prepare for output at this solution step.
     *
     * @param[in] t Time associated with field.
     * @param[in] mesh Mesh for output.
     */
    virtual
    void _openDataStep(const PylithReal t,
                       const pylith::topology::Mesh& mesh);

    /// Finalize output at this solution step.
    virtual
    void _closeDataStep(void);

    /** Write output for step in solution.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     */
    virtual
    void _writeDataStep(const PylithReal t,
                        const PylithInt tindex,
                        const pylith::topology::Field& solution);

    /** Append finite-element vertex field to file.
     *
     * @param[in] t Time associated with field.
     * @param[in] field Field to output.
     * @param[in] mesh Mesh for output.
     */
    virtual
    void _appendField(const PylithReal t,
                      pylith::topology::Field* field,
                      const pylith::topology::Mesh& mesh);

    /** Names of information fields for output.
     *
     * Expand "all" into list of actual fields.
     *
     * @param[in] auxField Auxiliary field.
     */
    pylith::string_vector _infoNamesExpanded(const pylith::topology::Field* auxField) const;

    /** Names of data fields for output.
     *
     * Expand "all" into list of actual fields in the solution field and derived field. We do not include subfields of
     * the auxiliary field, because most are physical parameters that do not change.
     *
     * @param[in] solution Solution field.
     * @param[in] auxField Auxiliary field.
     * @param[in] derivedField Derived field.
     */
    pylith::string_vector _dataNamesExpanded(const pylith::topology::Field& solution,
                                             const pylith::topology::Field* auxField,
                                             const pylith::topology::Field* derivedField) const;

    /** Get buffer for field.
     *
     * Find the most appropriate buffer that matches field, reusing and reallocating as necessary.
     *
     * @param[in] fieldIn Input field.
     * @param[in] name Name of subfield (optional).
     * @returns Field to use as buffer for outputting field.
     */
    pylith::topology::Field* _getBuffer(const pylith::topology::Field& fieldIn,
                                        const char* name=NULL);

    /** Dimension field.
     *
     * @param[in] fieldIn Field to dimensionalize.
     */
    pylith::topology::Field* _dimension(pylith::topology::Field* fieldIn);

    /** Get basis order of field.
     *
     * @param[in] field Field with one subfield for output.
     *
     * @returns Basis order if field contains single subfield, otherwise -1;
     */
    int _basisOrder(const pylith::topology::Field& field);

    /** TEMPOARY Set label and label id.
     *
     * @param[in] label Name of label.
     * @param[in] labelId Label value.
     */
    void _temporarySetLabel(const char* label,
                            const PylithInt labelId);

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    pylith::topology::Fields* _fields;   ///< Container with field buffers used for output.
    DataWriter* _writer;   ///< Writer for data.
    FieldFilter* _fieldFilter;   ///< Filter applied to fields.
    OutputTrigger* _trigger; ///< Trigger for deciding how often to write output.

    pylith::string_vector _infoFields;
    pylith::string_vector _dataFields;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    // :TODO: Remove _label and _labelId once materials use their own PetscDM.
    std::string _label; ///< Name of label defining cells to include in output (=0 means use all cells in mesh).
    PylithInt _labelId; ///< Value of label defining which cells to include.

    static const char* _pyreComponent; ///< Name of Pyre component.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    OutputManager(const OutputManager&);   ///< Not implemented.
    const OutputManager& operator=(const OutputManager&);   ///< Not implemented

}; // OutputManager

#endif // pylith_meshio_OutputManager_hh


// End of file
