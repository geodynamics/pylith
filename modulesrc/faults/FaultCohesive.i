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

/** @file modulesrc/faults/FaultCohesive.i
 *
 * @brief Python interface to C++ FaultCohesive object.
 */

namespace pylith {
    namespace faults {
	
	class FaultCohesive : public pylith::feassemble::IntegratorPointwise {

	    // PUBLIC METHODS /////////////////////////////////////////////////
	public :

	    /// Default constructor.
	    FaultCohesive(void);
      
	    /// Destructor.
	    virtual
	    ~FaultCohesive(void);
	    
	    /// Deallocate PETSc and local data structures.
	    virtual
	    void deallocate(void);
	    
	    /** Set material identifier of fault.
	     *
	     * @param[in] value Fault identifier
	     */
	    void id(const int value);
	    
	    /** Get material identifier of fault.
	     *
	     * @returns Fault identifier
	     */
	    int id(void) const;
	    
	    /** Set label of group of vertices associated with fault.
	     *
	     * @param[in] value Label of fault
	     */
	    void label(const char* value);
	    
	    /** Get label of group of vertices associated with fault.
	     *
	     * @returns Label of fault
	     */
	    const char* label(void) const;
	    
	    /** Set label of group of vertices defining buried edge of fault.
	     *
	     * @param[in] value Label of fault
	     */
	    void edge(const char* value);
	    
	    /** Get label of group of vertices defining buried edge of fault.
	     *
	     * @returns Label of fault
	     */
	    const char* edge(void) const;
	    
	    /** Set first choice for reference direction to discriminate among tangential directions in 3-D.
	     *
	     * @param vec Reference direction unit vector.
	     */
	    void refDir1(const PylithReal vec[3]);
	    
	    /** Set second choice for reference direction to discriminate among tangential directions in 3-D.
	     *
	     * @param vec Reference direction unit vector.
	     */
	    void refDir2(const PylithReal vec[3]);
	    
	    /** Adjust mesh topology for fault implementation.
	     *
	     * @param mesh[in] PETSc mesh.
	     */
	    void adjustTopology(pylith::topology::Mesh* const mesh);
	    
	    /** Verify configuration is acceptable.
	     *
	     * @param[in] solution Solution field.
	     */
	    virtual
	    void verifyConfiguration(const pylith::topology::Field& solution) const;
	    
	    /** Initialize fault.
	     *
	     * Create fault mesh from cohesive cells and cohesive point map.
	     *
	     * Derived class initialize, should:
	     * 1. Setup subfields in auxiliary field.
	     * 2. Populate auxiliary subfields.
	     * 3. Set finite-element kernels.
	     *
	     * @param[in] solution Solution field (layout).
	     */
	    virtual
	    void initialize(const pylith::topology::Field& solution);
	    
	    
	}; // class FaultCohesive
	
    } // faults
} // pylith


// End of file 
