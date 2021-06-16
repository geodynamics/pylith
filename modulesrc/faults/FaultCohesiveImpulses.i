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

/** @file modulesrc/faults/FaultCohesiveKin.i
 *
 * @brief Python interface to C++ FaultCohesiveKin object.
 */

namespace pylith {
    namespace faults {
        class pylith::faults::FaultCohesiveImpulses : public pylith::faults::FaultCohesiveKin {
            friend class TestFaultCohesiveImpulses; // unit testing

            // PUBLIC METHODS
            // //////////////////////////////////////////////////////////////////////////////////////////////////
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

            // PROTECTED METHODS
            // ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

            /** Update slip subfield in auxiliary field at beginning of time step.
             *
             * @param[out] auxiliaryField Auxiliary field.
             * @param[in] t Current time.
             */
            void _updateSlip(pylith::topology::Field* auxiliaryField,
                             const double t);

        }; // class FaultCohesiveImpulses

    } // faults
} // pylith

// End of file
