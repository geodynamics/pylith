#!/usr/bin/env python3
"""Application for uploading release assets to GitHub repository.
"""

import argparse
import requests
import mimetypes
import os

class App:
    BASE_URL = "https://api.github.com/repos/geodynamics/{package}/releases"

    def set_token(self):
        import getpass

        token = getpass.getpass(prompt="Enter GitHub personal access token: ")
        os.environ["GITHUB_TOKEN"] = token

    def upload_file(self,
        filename: str,
        release: str,
        package: str="pylith",
        ):

        self._get_token()
        self._get_release_info(package, release)

        mimetype = mimetypes.guess_type(filename)[0]
        if mimetype is None:
            raise IOError(f"Unknown mimetype for '{filename}'.")

        headers = {
            "Authorization": f"token {self.token}",
            "Content-Type": mimetype,
            "Accept": "application/vnd.github.v3+json",
            }
        params = (("name", filename),)
        data = open(filename, "rb").read()
        response = requests.post(self.upload_url, headers=headers, params=params, data=data)
        if response.status_code == 201:
            download_url = response.json()["browser_download_url"]
            print(f"Upload complete. Download URL: {download_url}")
        else:
            print("Error while uploading file.")
            print(response.text)

    def _get_token(self):
        if "GITHUB_TOKEN" not in os.environ:
            raise IOError("Could not find GitHub token in environment.")
        self.token = os.environ["GITHUB_TOKEN"]

    def _get_release_info(self, package, tag):
        url = self.BASE_URL.format(package=package)
        headers = {"Authorization": f"token {self.token}"}

        response = requests.get(url, headers=headers)
        releases = response.json()
        for release in releases:
            if release["tag_name"] == tag:
                self.release_id = release["id"]
                self.upload_url = release["upload_url"].replace("{?name,label}","")
                return
        tags = "\n\t".join([release["tag_name"] for release in releases])
        raise ValueError(
            f"Could not find release '{tag}' in {package} GitHub repository.\n"
            f"Releases:\n\t{tags}"
            )

def cli():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--set-token", dest="set_token", action="store_true", help="GitHub token")
    parser.add_argument("--package", dest="package", action="store", default="pylith", help="GitHub repository")
    parser.add_argument("--release", dest="release", action="store", required=True, help="Release tag")
    parser.add_argument("--asset", dest="asset", action="store", required=True, help="File to upload")
    args = parser.parse_args()

    app = App()
    if args.set_token:
        app.set_token()
    if args.asset:
        app.upload_file(filename=args.asset, release=args.release, package=args.package)

if __name__ == "__main__":
    cli()
