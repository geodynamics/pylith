// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/TimeDependentAuxiliaryFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary fields for time-dependent boundary conditions.
 */

#if !defined(pylith_bc_timedependentauxiliaryfactory_hh)
#define pylith_bc_timedependentauxiliaryfactory_hh

// Include directives ---------------------------------------------------
#include "bcfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

// TimeDependentAuxiliaryFactory-----------------------------------------------
/// @brief C++ helper class for setting up auxiliary fields for time-dependent boundary conditions.
class pylith::bc::TimeDependentAuxiliaryFactory : public pylith::feassemble::AuxiliaryFactory {
    friend class TestDirichletAuxiliaryFactory;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    TimeDependentAuxiliaryFactory(void);

    /// Destructor.
    ~TimeDependentAuxiliaryFactory(void);

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Add initial amplitude field to auxiliary fields.
    void initialAmplitude(void);

    /// Add rate amplitude field to auxiliary fields.
    void rateAmplitude(void);

    /// Add rate start time amplitude field to auxiliary fields.
    void rateStartTime(void);

    /// Add time history amplitude field to auxiliary fields.
    void timeHistoryAmplitude(void);

    /// Add time history start time field to auxiliary fields.
    void timeHistoryStartTime(void);

    /// Add time history value field to auxiliary fields.
    void timeHistoryValue(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    static const char* _genericComponent; ///< Name of generic component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    TimeDependentAuxiliaryFactory(const TimeDependentAuxiliaryFactory &);   ///< Not implemented.
    const TimeDependentAuxiliaryFactory& operator=(const TimeDependentAuxiliaryFactory&);   ///< Not implemented

}; // class TimeDependentAuxiliaryFactory

#endif // pylith_bc_timedependentauxiliaryfactory_hh


// End of file
