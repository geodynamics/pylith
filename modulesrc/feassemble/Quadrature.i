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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/** @file modulesrc/feassemble/Quadrature.i
 *
 * @brief Python interface to C++ Quadrature object.
 */

namespace pylith {
  namespace feassemble {

    template<typename mesh_type>
    class Quadrature : public QuadratureRefCell
    { // Quadrature

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      Quadrature(void);
      
      /// Destructor
      ~Quadrature(void);
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Copy constructor.
       *
       * @param q Quadrature to copy
       */
      Quadrature(const Quadrature& q);
      
      /** Set flag for checking ill-conditioning.
       *
       * @param flag True to check for ill-conditioning, false otherwise.
       */
      void checkConditioning(const bool flag);
      
      /** Get flag for checking ill-conditioning.
       *
       * @returns True if checking for ill-conditioning, false otherwise.
       */
      bool checkConditioning(void) const;

      /// Setup quadrature engine.
      void initializeGeometry(void);
      
      /// Deallocate temporary storage.
      void clear(void);

      /** Get precomputed geometry fields.
       *
       * @returns Geometry fields.
       */
      const pylith::topology::Fields<pylith::topology::Field<mesh_type> >* geometryFields(void) const;

    }; // Quadrature

  } // feassemble
} // pylith


// End of file 
