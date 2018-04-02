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
 * @file modulesrc/meshio/OutputSoln.i
 *
 * @brief Python interface to C++ OutputSoln object.
 */

namespace pylith {
    namespace meshio {

	class pylith::meshio::OutputSoln : public OutputManager {

	// PUBLIC METHODS ///////////////////////////////////////////////
	public :
	    
	    /// Constructor
	    OutputSoln(void);
	    
	    /// Destructor
	    ~OutputSoln(void);
	    
	    /// Deallocate PETSc and local data structures.
	    void deallocate(void);
	    
	    /** Set names of solution fields to output.
	     *
	     * @param[in] names Array of names of fields to output.
	     * @param[in] numNames Length of array.
	     */
	    %apply(const char* const* string_list, const int list_len){
			(const char* names[], const int numNames)
	    };
	    void vertexDataFields(const char* names[],
				  const int numNames);
        %clear(const char* const* names, const int numNames);

	    /** Verify configuration.
	     *
	     * @param[in] solution Solution field.
	     * @param[in] auxField Auxiliary field.
	     */
	    void verifyConfiguration(const pylith::topology::Field& solution,
				     const pylith::topology::Field& auxField) const;
	    
	    /** Write solution at time step.
	     *
	     * @param[in] t Current time.
	     * @param[in] tindex Current time step.
	     * @param[in] solution Solution at time t.
	     * @param[in] auxField Auxiliary field.
	     */
	    virtual
	    void writeTimeStep(const PylithReal t,
			       const PylithInt tindex,
			       const pylith::topology::Field& solution,
			       const pylith::topology::Field& auxField);
	    
	}; // OutputSoln
      
    } // meshio
} // pylith


// End of file
