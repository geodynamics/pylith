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
 * @file pylith/utils/EventLogger.hh
 *
 * @brief C++ object for managing event logging using PETSc.
 *
 * Each logger object manages the events for a single "logging class".
 */

#if !defined(pylith_utils_eventlogger_hh)
#define pylith_utils_eventlogger_hh

#include <string> // USES std::string
#include <map> // USES std::map

#include "petsclog.h" // USES PetscLogEventBegin/End() in inline methods

namespace pylith {
  namespace utils {
    class EventLogger;
    class TestEventLogger; // unit testing
  } // utils
} // pylith

class pylith::utils::EventLogger
{ // EventLogger
  friend class TestEventLogger; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  EventLogger(void);

  /// Destructor
  ~EventLogger(void);

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

  /// Setup logging class.
  void initialize(void);

  /** Register event.
   *
   * @prerequisite Must call initialize() before registerEvent().
   * 
   * @param name Name of event.
   * @returns Event identifier.
   */
  int registerEvent(const char* name);

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
  void eventBegin(const int id);

  /** Log event end.
   *
   * @param id Event identifier.
   */
  void eventEnd(const int id);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  EventLogger(const EventLogger&); ///< Not implemented
  const EventLogger& operator=(const EventLogger&); ///< Not implemented

// PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private :

  typedef std::map<std::string,int> map_event_type;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  std::string _className; ///< Name of logging class
  int _classId; ///< PETSc logging identifier for class
  map_event_type _events; ///< PETSc logging identifiers for events

}; // EventLogger

#include "EventLogger.icc" // inline methods

#endif // pylith_utils_eventlogger_hh


// End of file 
