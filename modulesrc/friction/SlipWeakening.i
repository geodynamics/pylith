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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/friction/SlipWeakening.i
 *
 * Python interface to C++ SlipWeakening object.
 */

namespace pylith {
  namespace friction {

    class SlipWeakening : public FrictionModel
    { // class SlipWeakening

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      SlipWeakening(void);

      /// Destructor.
      ~SlipWeakening(void);

      /** Compute properties from values in spatial database.
       *
       * @param flag True if forcing healing, false otherwise.
       */
      void forceHealing(const bool flag);

      // PROTECTED METHODS //////////////////////////////////////////////
    protected :

      /** Compute properties from values in spatial database.
       *
       * @param propValues Array of property values.
       * @param dbValues Array of database values.
       */
      void _dbToProperties(PylithScalar* const propValues,
			   const scalar_array& dbValues) const;

      /** Nondimensionalize properties.
       *
       * @param values Array of property values.
       * @param nvalues Number of values.
       */
      void _nondimProperties(PylithScalar* const values,
			     const int nvalues) const;

      /** Dimensionalize properties.
       *
       * @param values Array of property values.
       * @param nvalues Number of values.
       */
      void _dimProperties(PylithScalar* const values,
			  const int nvalues) const;

      /** Compute friction from properties and state variables.
       *
       * @param t Current time.
       * @param slip Current slip at location.
       * @param slipRate Current slip rate at location.
       * @param normalTraction Normal traction at location.
       * @param properties Properties at location.
       * @param numProperties Number of properties.
       * @param stateVars State variables at location.
       * @param numStateVars Number of state variables.
       */
      PylithScalar _calcFriction(const PylithScalar t,
				 const PylithScalar slip,
				 const PylithScalar slipRate,
				 const PylithScalar normalTraction,
				 const PylithScalar* properties,
				 const int numProperties,
				 const PylithScalar* stateVars,
				 const int numStateVars);

      /** Update state variables (for next time step).
       *
       * @param slip Current slip at location.
       * @param slipRate Current slip rate at location.
       * @param normalTraction Normal traction at location.
       * @param stateVars State variables at location.
       * @param numStateVars Number of state variables.
       * @param properties Properties at location.
       * @param numProperties Number of properties.
       */
      void _updateStateVars(const PylithScalar slip,
			    const PylithScalar slipRate,
			    const PylithScalar normalTraction,
			    PylithScalar* const stateVars,
			    const int numStateVars,
			    const PylithScalar* properties,
			    const int numProperties);
      
    }; // class SlipWeakening

  } // friction
} // pylith


// End of file 
