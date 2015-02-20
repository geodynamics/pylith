#!/usr/bin/env python

url = "https://api.github.com/repos/geodynamics/pylith/releases"

import urllib2
raw = urllib2.urlopen(url).read()

import json
data = json.loads(raw)
count = 0
for release in data:
    print release['name']
    for asset in release['assets']:
        print "    ",asset['name'],asset['download_count']
        if not "petsc" in asset['name'] and not "manual" in asset['name'] and not "installer" in asset['name']:
            count += asset['download_count']
print "Total downloads:",count
