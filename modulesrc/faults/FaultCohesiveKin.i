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
        class FaultCohesiveKin : public pylith::faults::FaultCohesive {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            FaultCohesiveKin(void);

            /// Destructor.
            virtual ~FaultCohesiveKin(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set kinematic earthquake sources.
             *
             * @param names Array of kinematic earthquake source names.
             * @param numNames Number of earthquake sources.
             * @param sources Array of kinematic earthquake sources.
             * @param numSources Number of earthquake sources.
             */
            %apply(const char* const* string_list, const int list_len) {
                (const char* const* names,
                 const int numNames)
            };
            void setEqRuptures(const char* const* names,
                               const int numNames,
                               pylith::faults::KinSrc** ruptures,
                               const int numRuptures);

            %clear(const char* const* names, const int numNames);

            /** Verify configuration is acceptable.
             *
             * @param[in] solution Solution field.
             */
            void verifyConfiguration(const pylith::topology::Field& solution) const;

            /** Create integrator and set kernels.
             *
             * @param[in] solution Solution field.
             * @param[in] materials Materials in problem.
             * @returns Integrator if applicable, otherwise NULL.
             */
            pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution,
                                                             const std::vector<pylith::materials::Material*>& materials);

            /** Create auxiliary field.
             *
             * @param[in] solution Solution field.
             * @param[in\ domainMesh Finite-element mesh associated with integration domain.
             *
             * @returns Auxiliary field if applicable, otherwise NULL.
             */
            pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                          const pylith::topology::Mesh& domainMesh);

            /** Update auxiliary subfields at beginning of time step.
             *
             * @param[out] auxiliaryField Auxiliary field.
             * @param[in] t Current time.
             */
            void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                      const double t);

            // PROTECTED METHODS //////////////////////////////////////////////////////////////////
protected:

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

            /** Update slip related subfields in auxiliary field at beginning of time step.
             *
             * @param[out] auxiliaryField Auxiliary field.
             * @param[in] t Current time.
             * @param[in] bitSlipSubfields Slip subfields to update.
             */
            void _updateSlip(pylith::topology::Field* auxiliaryField,
                             const double t,
                             const int bitSlipSubfields);

            /** Set kernels for residual.
             *
             * @param[out] integrator Integrator for material.
             * @param[in] solution Solution field.
             * @param[in] materials Materials in problem.
             */
            void _setKernelsResidual(pylith::feassemble::IntegratorInterface* integrator,
                                     const pylith::topology::Field& solution,
                                     const std::vector<pylith::materials::Material*>& materials) const;

            /** Set kernels for Jacobian.
             *
             * @param[out] integrator Integrator for material.
             * @param[in] solution Solution field.
             * @param[in] materials Materials in problem.
             */
            void _setKernelsJacobian(pylith::feassemble::IntegratorInterface* integrator,
                                     const pylith::topology::Field& solution,
                                     const std::vector<pylith::materials::Material*>& materials) const;

        }; // class FaultCohesiveKin

    } // faults
} // pylith

// End of file
