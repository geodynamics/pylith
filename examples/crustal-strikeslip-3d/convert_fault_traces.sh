#!/bin/bash

for fault in main east west; do
    cat faulttrace_${fault}_latlon.txt | cs2cs epsg:4326 epsg:32611 > faulttrace_${fault}_utm.txt
done
