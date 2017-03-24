#!/usr/bin/env python

packages = ["pylith","pylith_installer","spatialdata"]
baseurl = "https://api.github.com/repos/geodynamics/%s/releases"

import urllib2
import json

for package in packages:
    raw = urllib2.urlopen(baseurl % package).read()

    data = json.loads(raw)
    count = 0
    print("\n%s" % package)
    for release in data:
        print("    %s" % release['name'])
        for asset in release['assets']:
            print("        %s: %d" % (asset['name'],asset['download_count']))
            if package == "pylith":
                if not "petsc" in asset['name'] and not "manual" in asset['name'] and not "installer" in asset['name']:
                    count += asset['download_count']
            else:
                count += asset['download_count']
    print("  Total downloads: %d" % count)
