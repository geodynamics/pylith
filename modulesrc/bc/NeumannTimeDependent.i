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

/** @file modulesrc/bc/NeumannTimeDependent.i
 *
 * @brief Python interface to C++ NeumannTimeDependent object.
 */

namespace pylith {
    namespace bc {

        class NeumannTimeDependent : public NeumannNew
        { // class NeumannTimeDependent

        // PUBLIC METHODS /////////////////////////////////////////////////
public:

        /// Default constructor.
        NeumannTimeDependent(void);

        /// Destructor.
        ~NeumannTimeDependent(void);

        /// Deallocate PETSc and local data structures.
        void deallocate(void);

	/** Set time history database.
	 *
	 * @param[in] db Time history database.
	 */
	void dbTimeHistory(spatialdata::spatialdb::TimeHistory* th);
	
	/** Get time history database.
	 *
	 * @preturns Time history database.
	 */
	const spatialdata::spatialdb::TimeHistory* dbTimeHistory(void);
	
	/** Use initial value term in time history expression.
	 *
	 * @param[in] value True if using initial value term in expression.
	 */
	void useInitial(const bool value);
	
	/** Get flag associated with using initial value term in time history expression.
	 *
	 * @returns True if using initial value term in expression, false otherwise.
	 */
	bool useInitial(void) const;
	
	/** Use rate value term in time history expression.
	 *
	 * @param[in] value True if using rate value term in expression.
	 */
	void useRate(const bool value);
	
	/** Get flag associated with using rate value term in time history expression.
	 *
	 * @returns True if using rate value term in expression, false otherwise.
	 */
	bool useRate(void) const;
	
	/** Use time history term in time history expression.
	 *
	 * @param[in] value True if using time history term in expression.
	 */
	void useTimeHistory(const bool value);
	
	/** Get flag associated with using time history term in time history expression.
	 *
	 * @returns True if using time history term in expression, false otherwise.
	 */
	bool useTimeHistory(void) const;
	
	/** Update auxiliary fields at beginning of time step.
	 *
	 * @param[in] t Current time.
	 * @param[in] dt Current time step.
	 */
	void prestep(const double t,
		     const double dt);
	
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
	 * order of the auxiliary fields in the FE kernels.
	 *
	 * @param[in] solution Solution field.
	 */
	void _auxFieldSetup(const pylith::topology::Field& solution);
	
	/** Set kernels for RHS residual G(t,s).
	 *
	 * Potentially, there are g0 and g1 kernels for each equation. If no
	 * kernel is needed, then set the kernel function to NULL.
	 *
	 * @param solution Solution field.
	 */
	void _setFEKernelsRHSResidual(const pylith::topology::Field& solution) const;
	
	/** Get factory for setting up auxliary fields.
	 *
	 * @returns Factor for auxiliary fields.
	 */
	pylith:feassemble::AuxiliaryFactory* _auxFactory(void);


        }; // class NeumannTimeDependent

    } // bc
} // pylith


// End of file
