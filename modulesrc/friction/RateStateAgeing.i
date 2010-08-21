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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/friction/RateStateAgeing.i
 *
 * Python interface to C++ RateStateAgeing object.
 */

namespace pylith {
  namespace friction {

    class RateStateAgeing : public FrictionModel
    { // class RateStateAgeing

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      RateStateAgeing(void);

      /// Destructor.
      ~RateStateAgeing(void);

      // PROTECTED METHODS //////////////////////////////////////////////
    protected :

      /** Compute properties from values in spatial database.
       *
       * @param propValues Array of property values.
       * @param dbValues Array of database values.
       */
      void _dbToProperties(double* const propValues,
			   const double_array& dbValues) const;

      /** Nondimensionalize properties.
       *
       * @param values Array of property values.
       * @param nvalues Number of values.
       */
      void _nondimProperties(double* const values,
			     const int nvalues) const;

      /** Dimensionalize properties.
       *
       * @param values Array of property values.
       * @param nvalues Number of values.
       */
      void _dimProperties(double* const values,
			  const int nvalues) const;

      /** Compute friction from properties and state variables.
       *
       * @param slip Current slip at location.
       * @param slipRate Current slip rate at location.
       * @param normalTraction Normal traction at location.
       * @param properties Properties at location.
       * @param numProperties Number of properties.
       * @param stateVars State variables at location.
       * @param numStateVars Number of state variables.
       */
      double _calcFriction(const double slip,
			   const double slipRate,
			   const double normalTraction,
			   const double* properties,
			   const int numProperties,
			   const double* stateVars,
			   const int numStateVars);

    }; // class RateStateAgeing

  } // friction
} // pylith


// End of file 
