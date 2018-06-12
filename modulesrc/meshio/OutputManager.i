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
 * @file modulesrc/meshio/OutputManager.i
 *
 * @brief Python interface to C++ OutputManager object.
 */

namespace pylith {
    namespace meshio {

        class pylith::meshio::OutputManager : public pylith::utils::PyreComponent,
					      public pylith::feassemble::Observer {

	  // PUBLIC METHODS /////////////////////////////////////////////////
	public :
	  
	  /// Constructor
	  OutputManager(void);
	  
	  /// Destructor
	  virtual
	  ~OutputManager(void);
	  
	  /// Deallocate PETSc and local data structures.
	  virtual
	  void deallocate(void);
	  
	  /** Set output trigger for how often to write output.
	   *
	   * @param[in] otrigger Output trigger.
	   */
	  void trigger(pylith::meshio::OutputTrigger* const otrigger);
	  
	  /** Get trigger for how often to write otuput.
	   *
	   * @returns Output trigger.
	   */
	  const pylith::meshio::OutputTrigger* trigger(void) const;
	  
	  /** Set writer to write data to file.
	   *
	   * @param datawriter Writer for data.
	   */
	  void writer(DataWriter* const datawriter);
	  
	  /** Set filter for data.
	   *
	   * @param filter Filter to apply to data before writing.
	   */
	  void fieldFilter(FieldFilter* const filter);
	  
	    /** Set names of information fields to output.
	     *
	     * @param[in] names Array of field names.
	     * @param[in] numNames Length of array.
	     */
	    %apply(const char* const* string_list, const int list_len){
			(const char* names[], const int numNames)
	    };
	    void infoFields(const char* names[],
			    const int numNames);
	    %clear(const char* const* names, const int numNames);

	    /** Set names of data fields to output.
	     *
	     * @param[in] names Array of field names.
	     * @param[in] numNames Length of array.
	     */
	    %apply(const char* const* string_list, const int list_len){
			(const char* names[], const int numNames)
	    };
	    void dataFields(const char* names[],
			    const int numNames);
	    %clear(const char* const* names, const int numNames);

	    /** Verify configuration.
	     *
	     * @param[in] solution Solution field.
	     */
	    virtual
	    void verifyConfiguration(const pylith::topology::Field& solution) const = 0;

	      /** Receive update from subject.
	       *
	       * @param[in] t Current time.
	       * @param[in] tindex Current time step.
	       * @param[in] solution Solution at time t.
	       * @param[in] infoOnly Flag is true if this update is before solution is available (e.g., after initialization).
	       */
	  virtual
	  void update(const PylithReal t,
		      const PylithInt tindex,
		      const pylith::topology::Field& solution,
		      const bool infoOnly=false);

	}; // OutputManager

    } // meshio
} // pylith


// End of file
