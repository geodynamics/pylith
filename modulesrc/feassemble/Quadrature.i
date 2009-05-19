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

      /// Setup quadrature engine.
      void initializeGeometry(void);
      
      /// Deallocate temporary storage.
      void clear(void);

      /** Get precomputed coordinates of quadrature points
       *
       * @returns Array of coordinates of quadrature points in cell
       */
      const pylith::topology::Field<mesh_type>& quadPtsPrecomp(void) const;

      /** Get precomputed derivatives of basis fns evaluated at quadrature points.
       *
       * @returns Array of derivatives of basis fns evaluated at
       * quadrature points
       */
      const pylith::topology::Field<mesh_type>& basisDerivPrecomp(void) const;

      /** Get precomputed Jacobians evaluated at quadrature points.
       *
       * @returns Array of Jacobian inverses evaluated at quadrature points.
       */
      const pylith::topology::Field<mesh_type>& jacobianPrecomp(void) const;

      /** Get precomputed determinants of Jacobian evaluated at quadrature points.
       *
       * @returns Array of determinants of Jacobian evaluated at quadrature pts
       */
      const pylith::topology::Field<mesh_type>& jacobianDetPrecomp(void) const;
      
    }; // Quadrature

  } // feassemble
} // pylith


// End of file 
