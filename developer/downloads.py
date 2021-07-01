#!/usr/bin/env python3

import json
import ssl
import urllib.request


packages = ["pylith", "pylith_installer", "spatialdata"]
baseurl = "https://api.github.com/repos/geodynamics/%s/releases"

ssl._create_default_https_context = ssl._create_unverified_context

for package in packages:
    raw = urllib.request.urlopen(baseurl % package).read()

    data = json.loads(raw)
    count = 0
    print("\n%s" % package)
    for release in data:
        print("    %s" % release['name'])
        for asset in release['assets']:
            print("        %s: %d" % (asset['name'], asset['download_count']))
            if package == "pylith":
                if not "petsc" in asset['name'] and not "manual" in asset['name'] and not "installer" in asset['name']:
                    count += asset['download_count']
            else:
                count += asset['download_count']
    print("  Total downloads: %d" % count)
