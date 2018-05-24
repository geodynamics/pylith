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

/** @file libsrc/faults/KinSrcAuxiliaryFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary fields for slip time functions.
 */

#if !defined(pylith_faults_kinsrcauxiliaryfactory_hh)
#define pylith_faults_kinsrcauxiliaryfactory_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

// KinSrcAuxiliaryFactory-----------------------------------------------
/// @brief C++ helper class for setting up auxiliary fields for faults.
class pylith::faults::KinSrcAuxiliaryFactory : public pylith::feassemble::AuxiliaryFactory {
    friend class TestSlipFnAuxiliaryFactory;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    KinSrcAuxiliaryFactory(void);

    /// Destructor.
    ~KinSrcAuxiliaryFactory(void);

    /// Add slip initiation time (relative to origin time) subfield to auxiliary fields.
    void initiationTime(void);

    /// Add rise time subfield to auxiliary fields.
    void riseTime(void);

    /// Add finalSlip to auxiliary fields.
    void finalSlip(void);

    /// Add slipRate to auxiliary fields.
    void slipRate(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    static const char* _genericComponent; ///< Name of generic component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    KinSrcAuxiliaryFactory(const KinSrcAuxiliaryFactory &);   ///< Not implemented.
    const KinSrcAuxiliaryFactory& operator=(const KinSrcAuxiliaryFactory&);   ///< Not implemented

}; // class KinSrcAuxiliaryFactory

#endif // pylith_faults_kinsrcauxiliaryfactory_hh


// End of file
