// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file modulesrc/meshio/OutputIntegrator.i
 *
 * @brief Python interface to C++ OutputIntegrator object.
 */

namespace pylith {
    namespace meshio {

	class pylith::meshio::OutputIntegrator : public OutputManager {

	// PUBLIC METHODS ///////////////////////////////////////////////
	public :
	    
	    /// Constructor
	    OutputIntegrator(void);
	    
	    /// Destructor
	    ~OutputIntegrator(void);
	    
	    /// Deallocate PETSc and local data structures.
	    void deallocate(void);
	    
	    /** Set names of information fields to output.
	     *
	     * @param[in] names Array of names of fields to output.
	     * @param[in] numNames Length of array.
	     */
	    %apply(const char* const* string_list, const int list_len){
			(const char* names[], const int numNames)
	    };
	    void vertexInfoFields(const char* names[],
				  const int numNames);
	    %clear(const char* const* names, const int numNames);

	    /** Verify configuration.
	     *
	     * @param[in] auxField Auxiliary field.
	     */
	    virtual
	    void verifyConfiguration(const pylith::topology::Field& auxField) const;
	    
	    /** Write information.
	     *
	     * @param[in] auxField Auxiliary field.
	     */
	    virtual
	    void writeInfo(const pylith::topology::Field& auxField);

	    /** Write solution at time step.
	     *
	     * @param[in] t Current time.
	     * @param[in] timeStep Current time step.
	     * @param[in] solution Solution at time t.
	     * @param[in] auxField Auxiliary field.
	     */
	  virtual
	    void writeTimeStep(const PylithReal t,
			       const PylithInt timeStep,
			       const pylith::topology::Field& solution,
			       const pylith::topology::Field& auxField);
	    
	}; // OutputIntegrator
      
    } // meshio
} // pylith


// End of file
