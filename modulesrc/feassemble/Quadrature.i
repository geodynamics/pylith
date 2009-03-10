// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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
      
      /// Deallocate temporary storage.
      void clear(void);
      
    }; // Quadrature

  } // feassemble
} // pylith


// End of file 
