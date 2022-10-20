#!/usr/bin/env python3

import requests

packages = ["pylith", "pylith_installer", "spatialdata"]
baseurl = "https://api.github.com/repos/geodynamics/{package}/releases"

for package in packages:
    url = baseurl.format(package=package)
    response = requests.get(url)
    data = response.json()
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
