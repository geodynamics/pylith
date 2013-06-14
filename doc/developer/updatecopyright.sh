#!/bin/bash

for suffix in .cc .icc .hh .h .am .py .i .jou; do
  for f in `find . -name "*$suffix"` configure.ac aclocal.m4 COPYING; do
    sed -e "s/Copyright (c) 2010-2011 University of California, Davis/Copyright (c) 2010-2013 University of California, Davis/g" $f > tmp && mv -f tmp $f
    sed -e "s/Copyright (c) 2010-2012 University of California, Davis/Copyright (c) 2010-2013 University of California, Davis/g" $f > tmp && mv -f tmp $f
  done
done
