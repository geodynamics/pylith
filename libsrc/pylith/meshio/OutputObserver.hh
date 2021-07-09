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
#include <map> // HASA std::map

class pylith::meshio::OutputObserver : public pylith::utils::PyreComponent {
    friend class TestOutputObserver; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
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

    /** Set basis order for output.
     *
     * @param[in] value Basis order for output.
     */
    void setOutputBasisOrder(const int value);

    /** Set time scale.
     *
     * @param[in] value Time scale for dimensionalizing time.
     */
    virtual
    void setTimeScale(const PylithReal value);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Set context.
     *
     * @param[in] mesh Mesh associated with output.
     */
    void _setContext(const pylith::topology::Mesh & mesh);

    /** Get output subfield, creating if necessary.
     *
     * @param[in] field Field containing subfields.
     * @param[in] submesh Submesh associated with output.
     * @param[in] name Name of subfield.
     */
    virtual
    OutputSubfield* _getSubfield(const pylith::topology::Field& field,
                                 const pylith::topology::Mesh& submesh,
                                 const char* name);

    /** Append subfield at current time to output.
     *
     * @param[in] t Current time.
     * @param[in] subfield Subfield to write.
     */
    void _appendField(const PylithReal t,
                      const pylith::meshio::OutputSubfield& subfield);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    PylithReal _timeScale; ///< Time scale for dimentionalizing time.
    std::map<std::string, OutputSubfield*> _subfields; ///< Subfields extracted for output.
    DataWriter* _writer; ///< Writer for data.
    OutputTrigger* _trigger; ///< Trigger for deciding how often to write output.
    int _outputBasisOrder; ///< Basis order for output.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    OutputObserver(const OutputObserver&); ///< Not implemented.
    const OutputObserver& operator=(const OutputObserver&); ///< Not implemented

}; // OutputObserver

#endif // pylith_meshio_outputobserver_hh

// End of file
