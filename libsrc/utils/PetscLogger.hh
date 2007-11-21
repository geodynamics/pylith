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

/**
 * @file pylith/utils/PetscLogger.hh
 *
 * @brief C++ object for managing PETSc logger.
 *
 * Each logger object manages the events for a single "logging class".
 */

#if !defined(pylith_utils_petsclogger_hh)
#define pylith_utils_petsclogger_hh

#include <string> // USES std::string

namespace pylith {
  namespace utils {
    class PetscLogger;
    class TestPetscLogger;
  } // utils
} // pylith

class pylith::utils::PetscLogger
{ // Integrator
  friend class TestPetscLogger; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  PetscLogger(void);

  /// Destructor
  ~PetscLogger(void);

  /** Set name of logging class.
   *
   * @param name Name of logging class.
   */
  void className(const char* name);

  /** Get name of logging class.
   *
   * @returns Name of logging class.
   */
  const char* className(void) const;

  /** Register event.
   *
   * @param name Name of event.
   */
  void registerEvent(const char* name);

  /** Get event identifier.
   *
   * @param name Name of event.
   * @returns Event identifier.
   */
  int eventId(const char* name);

  /** Log event begin.
   *
   * @param id Event identifier.
   */
  void eventBegin(const int id)

  /** Log event end.
   *
   * @param id Event identifier.
   */
  void eventEnd(const int id)


// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  PetscLogger(const PetscLogger&); ///< Not implemented
  const PetscLogger& operator=(const PetscLogger&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  std::string _className; ///< Name of logging class
  std::map<string, int> _events;

}; // PetscLogger

#include "PetscLogger.icc" // inline methods

#endif // pylith_utils_petsclogger_hh


// End of file 
