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

#include <portinfo>

#include "LineParser.hh" // implementation of class methods

#include <iostream>
#include <sstream>


// ----------------------------------------------------------------------
pylith::utils::LineParser::LineParser(std::istream& sin,
				      const char* delimiter,
				      const size_t bufsize) :
  _in(sin),
  _delimiter(delimiter),
  _bufsize(bufsize),
  _buffer(0)
{ // constructor
  _buffer = new char[bufsize];
} // constructor

// ----------------------------------------------------------------------
pylith::utils::LineParser::~LineParser(void)
{ // destructor
  delete[] _buffer; _buffer = 0;
} // destructor

// ----------------------------------------------------------------------
const std::string&
pylith::utils::LineParser::next(void) {
  _in.getline(_buffer, _bufsize);
  _value = _buffer;
  size_t pos = _value.find(_delimiter);
  while (0 == pos && !_in.eof() && _in.good()) {
    _in.getline(_buffer, _bufsize);
    _value = _buffer;
    pos = _value.find(_delimiter);
  } // while
  const size_t last = (pos == std::string::npos) ? _value.length() : pos;
  _value = _value.substr(0, last);
  return _value;
} // next

// ----------------------------------------------------------------------
// Ignore input until character read.
void
pylith::utils::LineParser::ignore(const char marker)
{ // ignore
  _in.ignore(_bufsize, marker);
} // ignore

// ----------------------------------------------------------------------
// Eat whitespace.
void
pylith::utils::LineParser::eatws(void)
{ // eatws
  _in >> std::ws;
} // eatws


// End of file 
