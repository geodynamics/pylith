#!/bin/bash

for suffix in .cc .icc .hh .c .h .am .py .i .jou .odb .ac .in .dat; do
  for f in `find . -name "*$suffix"` aclocal.m4 COPYING ; do
    sed -e "s/Copyright (c) 2010-2015 University of California, Davis/Copyright (c) 2010-2017 University of California, Davis/g" $f > tmp && mv -f tmp $f
    sed -e "s/Copyright (c) 2010-2016 University of California, Davis/Copyright (c) 2010-2017 University of California, Davis/g" $f > tmp && mv -f tmp $f
    sed -e "s/Copyright (c) 2010-2017 University of California, Davis/Copyright (c) 2010-2017 University of California, Davis/g" $f > tmp && mv -f tmp $f
    sed -e "s/Matthew G. Knepley, University of Chicago/Matthew G. Knepley, Rice University/g" $f > tmp && mv -f tmp $f
  done
done
