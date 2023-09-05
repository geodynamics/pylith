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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/faults/FaultCohesiveKin.i
 *
 * @brief Python interface to C++ FaultCohesiveKin object.
 */

namespace pylith {
    namespace faults {
        class pylith::faults::FaultCohesiveImpulses: public pylith::faults::FaultCohesiveKin {
            friend class TestFaultCohesiveImpulses; // unit testing

            // PUBLIC METHODS /////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            FaultCohesiveImpulses(void);

            /// Destructor.
            ~FaultCohesiveImpulses(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set indices of fault degrees of freedom associated with
             * impulses.
             *
             * @param flags Array of indices for degrees of freedom.
             * @param size Size of array
             */
            %apply(int* INPLACE_ARRAY1, int DIM1) {
                (const int* flags,
                 const size_t size)
            };
            void setImpulseDOF(const int* flags,
                               const size_t size);

            %clear(const int* flags, const size_t size);

            /** Set threshold for nonzero impulse amplitude.
             *
             * @param value Threshold for detecting nonzero amplitude.
             */
            void setThreshold(const double value);

            /** Get the total number of impulses that will be applied.
             *
             * @returns Number of impulses.
             */
            size_t getNumImpulsesLocal(void);

            /** Verify configuration is acceptable.
             *
             * @param[in] solution Solution field.
             */
            void verifyConfiguration(const pylith::topology::Field& solution) const;

            /** Create auxiliary field.
             *
             * @param[in] solution Solution field.
             * @param[in] domainMesh Finite-element mesh associated with integration domain.
             *
             * @returns Auxiliary field if applicable, otherwise NULL.
             */
            pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                          const pylith::topology::Mesh& domainMesh);

            /** Create derived field.
             *
             * @param[in] solution Solution field.
             * @param[in\ domainMesh Finite-element mesh associated with integration domain.
             *
             * @returns Derived field if applicable, otherwise NULL.
             */
            pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                        const pylith::topology::Mesh& domainMesh);

            // PROTECTED METHODS //////////////////////////////////////////////////////////////////
protected:

            /** Update slip subfield in auxiliary field at beginning of time step.
             *
             * @param[out] auxiliaryField Auxiliary field.
             * @param[in] impulseIndex Index of impulse.
             */
            void _updateSlip(pylith::topology::Field* auxiliaryField,
                             const long impulseIndex);

            /** Set kernels for residual.
             *
             * @param[out] integrator Integrator for material.
             * @param[in] solution Solution field.
             * @param[in] materials Materials in problem.
             */
            void _setKernelsResidual(pylith::feassemble::IntegratorInterface* integrator,
                                     const pylith::topology::Field& solution,
                                     const std::vector < pylith::materials::Material* >& materials) const;

            /** Set kernels for Jacobian.
             *
             * @param[out] integrator Integrator for material.
             * @param[in] solution Solution field.
             * @param[in] materials Materials in problem.
             */
            void _setKernelsJacobian(pylith::feassemble::IntegratorInterface* integrator,
                                     const pylith::topology::Field& solution,
                                     const std::vector < pylith::materials::Material* >& materials) const;

        }; // class FaultCohesiveImpulses

    } // faults
} // pylith

// End of file
