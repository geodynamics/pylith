// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/** @file modulesrc/sources/RickerWavelet.i
 *
 * Python interface to C++ RickerWavelet.
 */

namespace pylith {
    namespace sources {
        class RickerWavelet : public pylith::sources::SourceTimeFunctionMomentTensorForce {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            RickerWavelet(void);

            /// Destructor.
            ~RickerWavelet(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            pylith::sources::AuxiliaryFactoryMomentTensorForce* getAuxiliaryFactory(void);

            /** Add source time function subfields to auxiliary field.
             *
             * @param[inout] auxiliaryField Auxiliary field.
             */
            void addAuxiliarySubfields(void);

            /** Get g1v kernel for residual, G(t,s).
            *
            * @param[in] coordsys Coordinate system.
            *
            * @return residual kernel for g1v
            .
            */
            PetscPointFunc getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const;

        };

        // class RickerWavelet

    } // sources
} // pylith

// End of file
