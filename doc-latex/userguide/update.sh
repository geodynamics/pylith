#!/bin/bash

texfiles=`find . -name "*.tex"`
for i in $texfiles; do
    sed -f input.sed $i > $i.tmp && mv -f $i.tmp $i
done
