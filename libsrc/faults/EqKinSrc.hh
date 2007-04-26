// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/EqKinSrc.hh
 *
 * @brief C++ object for managing parameters for a kinematic
 * earthquake source.
 *
 * EqKinSrc is responsible for providing the value of slip at time t
 * over a fault surface.
 */

#if !defined(pylith_faults_eqkinsrc_hh)
#define pylith_faults_eqkinsrc_hh

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class EqKinSrc;
    class TestEqKinSrc; // unit testing

    class SlipTimeFn; // HOLDSA SlipTimeFn
  } // faults
} // pylith

/// Namespace for spatialdata package
namespace spatialdata {
  namespace spatialdb {
    class SpatialDB;
  } // spatialdb

  namespace geocoords {
    class CoordSys;
  } // geocoords
} // spatialdata

/// C++ oject for managing parameters for a kinematic earthquake source.
class pylith::faults::EqKinSrc
{ // class EqKinSrc
  friend class TestEqKinSrc; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  EqKinSrc(void);

  /// Destructor.
  virtual
  ~EqKinSrc(void);

  /** Create copy of fault.
   *
   * @returns Copy of fault.
   */
  virtual
  EqKinSrc* clone(void) const;

  /** Set slip time function.
   *
   * @param slipfn Slip time function.
   */
  void slipfn(SlipTimeFn* slipfn);

  /** Initialize slip time function.
   *
   * @param mesh Finite-element mesh.
   * @param faultMesh Finite-element mesh of fault.
   * @param vertices Vertices where slip will be prescribed.
   * @param cs Coordinate system for mesh
   */
  virtual
  void initialize(const ALE::Obj<Mesh>& mesh,
		  const ALE::Obj<Mesh>& faultMesh,
		  const std::vector<Mesh::point_type>& vertices,
		  const spatialdata::geocoords::CoordSys* cs);

  /** Get slip on fault surface at time t.
   *
   * @param t Time t.
   * @param vertices Vertices where slip will be prescribed.
   */
  virtual
  const ALE::Obj<real_section_type>& slip(const double t,
			      const std::vector<Mesh::point_type>& vertices);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param s Source to copy
   */
  EqKinSrc(const EqKinSrc& s);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const EqKinSrc& operator=(const EqKinSrc& s);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  SlipTimeFn* _slipfn; ///< Slip time function

}; // class EqKinSrc

#include "EqKinSrc.icc" // inline methods

#endif // pylith_faults_eqkinsrc_hh


// End of file 
