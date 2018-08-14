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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/bc/ConstraintBoundary.i
 *
 * @brief Python interface to C++ abstract ConstraintBoundary object.
 */

namespace pylith {
    namespace bc {

        class ConstraintBoundary : public pylith::feassemble::ConstraintPointwise {

            // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

            /// Default constructor.
            ConstraintBoundary(void);

            /// Destructor.
            ~ConstraintBoundary(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set mesh label associated with boundary condition surface.
             *
             * @param[in] value Label of surface (from mesh generator).
             */
            void label(const char* value);

            /** Get mesh label associated with boundary condition surface.
             *
             * @returns Label of surface (from mesh generator).
             */
            const char* label(void) const;

            /** Set name of field in solution to constrain.
             *
             * @param[in] value Name of field in solution to constrain.
             */
            void field(const char* value);

            /** Get name of field in solution to constrain.
             *
             * @returns Name of field in solution to constrain.
             */
            const char* field(void) const;

            /** Get mesh associated with integrator domain.
             *
             * @returns Mesh associated with integrator domain.
             */
            const pylith::topology::Mesh& domainMesh(void) const;

            /** Verify configuration is acceptable.
             *
             * @param[in] solution Solution field.
             */
            virtual
            void verifyConfiguration(const pylith::topology::Field& solution) const;

            /** Initialize boundary condition.
             *
             * @param[in] solution Solution field.
             */
            void initialize(const pylith::topology::Field& solution);

            /** Set constrained values in solution field.
             *
             * @param[out] solution Solution field.
             * @param[in] t Current time.
             */
            void setSolution(pylith::topology::Field* solution,
                             const double t);


            // PROTECTED METHODS //////////////////////////////////////////////////
protected:

            /** Setup auxiliary subfields (discretization and query fns).
             *
             * Create subfields in auxiliary fields (includes name of the field,
             * vector field type, discretization, and scale for
             * nondimensionalization) and set query functions for filling them
             * from a spatial database.
             *
             * @attention The order of the calls to subfieldAdd() must match the
             * order of the auxiliary subfields in the FE kernels.
             *
             * @param[in] solution Solution field.
             */
            virtual
            void _auxFieldSetup(const pylith::topology::Field& solution) = 0;

            /** Set kernel for computing value of constrained degree of freedom from auxiliary field.
             *
             * @param solution Solution field.
             */
            virtual
            void _setFEKernelConstraint(const pylith::topology::Field& solution) = 0;

            /** Get point-wise function (kernel) for settings constraint from auxiliary field.
             *
             * @returns Point-wise function.
             */
            virtual
            PetscPointFunc _getFEKernelConstraint(void) = 0;

        }; // ConstraintBoundary

    } // bc
} // pylith


// End of file
