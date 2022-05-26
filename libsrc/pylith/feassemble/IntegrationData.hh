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
 * @file libsrc/problems/IntegrationData.hh
 *
 * @brief Object for managing data needed to integrate governing equations.
 */

#if !defined(pylith_problems_integrationdata_hh)
#define pylith_problems_integrationdata_hh

#include "pylith/problems/problemsfwd.hh" // Forward declarations
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/topologyfwd.hh" // USES Field

#include <string> // USES std::string
#include <map> // USES std::map

class pylith::problems::IntegrationData : pylith::utils::GenericComponent {
    friend class TestIntegrationData; // unit testing

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    static const std::string time;
    static const std::string time_step;
    static const std::string s_tshift;
    static const std::string t_state;
    static const std::string dt_residual;
    static const std::string dt_jacobian;
    static const std::string dt_lumped_jacobian_inverse;

    static const std::string solution;
    static const std::string solution_dot;
    static const std::string residual;
    static const std::string lumped_jacobian_inverse;

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    IntegrationData(void);

    /// Destructor
    ~IntegrationData(void);

    /// Deallocate data.
    void deallocate(void);

    /** Set scalar quantity.
     *
     * @param[in] name Name of scalar.
     * @param[in] value Value of scalar.
     */
    void setScalar(const std::string& name,
                   const PylithReal value);

    /** Get scalar quantity.
     *
     * @param[in] name Name of scalar.
     * @returns Value of scalar.
     */
    PylithReal getScalar(const std::string& name) const;

    /** Remove scalar quantity.
     *
     * @param[in] name Name of scalar.
     */
    void removeScalar(const std::string& name);

    /** Check if we have field with given name.
     *
     * @param[in] name Name of field.
     * @returns True if we have field, otherwise false.
     */
    bool hasField(const std::string& name) const;

    /** Set field.
     *
     * @param[in] name Name of field.
     * @param[in] field Field over the domain.
     */
    void setField(const std::string& name,
                  pylith::topology::Field* const field);

    /** Get field.
     *
     * @param[in] name Name of field.
     * @returns Field over the domain.
     */
    pylith::topology::Field* getField(const std::string& name) const;

    /** Dump integration data to std::string.
     *
     * @returns Integration data information as string.
     */
    std::string str(void) const;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    typedef std::map<std::string, PylithReal> scalars_map_t;
    typedef std::map<std::string, pylith::topology::Field*> fields_map_t;

    scalars_map_t _scalars;
    fields_map_t _fields;

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    IntegrationData(const IntegrationData&); ///< Not implemented.
    const IntegrationData& operator=(const IntegrationData&); ///< Not implemented.

}; // IntegrationData

#endif // pylith_problems_integrationdata_hh

// End of file
