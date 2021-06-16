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

/** @file libsrc/faults/FaultCohesiveImpulses.hh
 *
 * @brief C++ implementation for a fault surface with kinematic
 * (prescribed) slip implemented with cohesive elements which we decompose into impulses for a Green function.
 */

#if !defined(pylith_faults_faultcohesiveimpulses_hh)
#define pylith_faults_faultcohesiveimpulses_hh

#include "FaultCohesiveKin.hh" // ISA FaultCohesiveKin

class pylith::faults::FaultCohesiveImpulses : public pylith::faults::FaultCohesiveKin {
    friend class TestFaultCohesiveImpulses; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    FaultCohesiveImpulses(void);

    /// Destructor.
    ~FaultCohesiveImpulses(void);

    /** Return the number of impulses
     *
     * @param[out] the number of impulses
     */
    PylithInt getNumImpulses(void);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Update slip subfield in auxiliary field at beginning of time step.
     *
     * @param[out] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     */
    void _updateSlip(pylith::topology::Field* auxiliaryField,
                     const double t);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    FaultCohesiveImpulses(const FaultCohesiveImpulses&); ///< Not implemented
    const FaultCohesiveImpulses& operator=(const FaultCohesiveImpulses&); ///< Not implemented.

}; // class FaultCohesiveImpulses

#endif // pylith_faults_faultcohesiveimpulses_hh

// End of file
