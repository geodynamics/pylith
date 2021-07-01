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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

// SWIG interface to C++ ViscousFriction object.

/* This is nearly identical to the C++ ViscousFriction header
 * file. There are a few important differences required by SWIG:
 *
 * (1) Instead of forward declaring the ViscousFriction class, we
 * embed the class definition within the namespace declarations.
 *
 * (2) We only include public members and methods and implementations
 * of abstract methods, because this is an interface file.
 */

namespace contrib {
  namespace friction {

    class ViscousFriction : public pylith::friction::FrictionModel
    { // class ViscousFriction

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      ViscousFriction(void);

      /// Destructor.
      ~ViscousFriction(void);

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
       * @param t Time in simulation.
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

      /** Compute derivative of friction with slip from properties and
       * state variables.
       *
       * @param t Time in simulation.
       * @param slip Current slip at location.
       * @param slipRate Current slip rate at location.
       * @param normalTraction Normal traction at location.
       * @param properties Properties at location.
       * @param numProperties Number of properties.
       * @param stateVars State variables at location.
       * @param numStateVars Number of state variables.
       *
       * @returns Derivative of friction (magnitude of shear traction) at vertex.
       */
      PylithScalar _calcFrictionDeriv(const PylithScalar t,
				      const PylithScalar slip,
				      const PylithScalar slipRate,
				      const PylithScalar normalTraction,
				      const PylithScalar* properties,
				      const int numProperties,
				      const PylithScalar* stateVars,
				      const int numStateVars);
    }; // class ViscousFriction

  } // friction
} // pylith


// End of file 
