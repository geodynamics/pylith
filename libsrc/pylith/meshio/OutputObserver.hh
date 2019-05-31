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
 * @file libsrc/meshio/OutputObserver.hh
 *
 * @brief Manager for output via observer.
 */

#if !defined(pylith_meshio_outputobserver_hh)
#define pylith_meshio_outputobserver_hh

#include "meshiofwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/array.hh" // HASA string_vector

class pylith::meshio::OutputObserver : public pylith::utils::PyreComponent {
    friend class TestOutputObserver; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    OutputObserver(void);

    /// Destructor
    virtual ~OutputObserver(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set output trigger for how often to write output.
     *
     * @param[in] trigger Output trigger.
     */
    void setTrigger(pylith::meshio::OutputTrigger* const trigger);

    /** Get trigger for how often to write otuput.
     *
     * @returns Output trigger.
     */
    const pylith::meshio::OutputTrigger* getTrigger(void) const;

    /** Set writer to write data to file.
     *
     * @param[in] datawriter Writer for data.
     */
    void setWriter(pylith::meshio::DataWriter* const writer);

    /** Set filter for vertex data.
     *
     * @param[in] filter Filter to apply to vertex data before writing.
     */
    void setFieldFilter(pylith::meshio::FieldFilter* const filter);

    /** Set time scale.
     *
     * @param[in] value Time scale for dimensionalizing time.
     */
    virtual
    void setTimeScale(const PylithReal value);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

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
    pylith::topology::Field* _dimensionField(pylith::topology::Field* fieldIn);

    /** Get basis order of field.
     *
     * @param[in] field Field with one subfield for output.
     *
     * @returns Basis order if field contains single subfield, otherwise -1;
     */
    int _getBasisOrder(const pylith::topology::Field& field);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    PylithReal _timeScale; ///< Time scale for dimentionalizing time.
    pylith::topology::Fields* _fields; ///< Container with field buffers used for output.
    DataWriter* _writer; ///< Writer for data.
    FieldFilter* _fieldFilter; ///< Filter applied to fields.
    OutputTrigger* _trigger; ///< Trigger for deciding how often to write output.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    OutputObserver(const OutputObserver&); ///< Not implemented.
    const OutputObserver& operator=(const OutputObserver&); ///< Not implemented

}; // OutputObserver

#endif // pylith_meshio_outputobserver_hh

// End of file
