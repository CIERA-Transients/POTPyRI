import requests
import re
import os
from astropy.time import Time
import astropy.utils.data

def make_dummy_wcs():
    header = {'WCSAXES': 2,
              'CTYPE1': 'RA---TAN',
              'CTYPE2': 'DEC--TAN',
              'CRPIX1': 1.0,
              'CRPIX2': 1.0,
              'CRVAL1': 30.0,
              'CRVAL2': -30.0,
              'CD1_1': -5.61222E-05,
              'CD1_2': 0.0,
              'CD2_1': 0.0,
              'CD2_2': 5.61222E-05}

    return(header)

def download_github_file(
        git_file_path: str,
        output_dir: str | None = None,
        base_url: str | None = None,
        use_cached : bool = True
    ) -> str | None:
    """
    See: https://gist.github.com/fkraeutli/66fa741d9a8c2a6a238a01d17ed0edc5#retrieving-lfs-files
    
    Args:
        git_file_path (str): Path to the file in the POTPyRI_test repository.
        output_dir (str, optional): Local directory to save the file (default: current directory).
        base_url (str, optional): Base URL of the POTPyRI_test repository
            Defaults to: "https://raw.githubusercontent.com/oirlab/HISPEC_Test_Data/main/"
        use_cached (bool, optional): Use cached file if it exists (default: True).
    
    Returns:
        str | None: Path to the downloaded file or None if download failed.
    """

    if base_url is None:
        base_url = "https://raw.githubusercontent.com/CIERA-Transients/POTPyRI_test/main/"
    
    # Output dir
    if output_dir is None:
        output_dir = astropy.utils.data._get_download_cache_loc()
    
    # Set the output path
    output_path = os.path.join(output_dir, os.path.basename(git_file_path))

    if os.path.exists(output_path) and use_cached:
        return output_path
    
    # Ensure the output directory exists
    if output_dir:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    try:
        # Get the file URL from GitHub's raw content
        url = base_url + git_file_path
        response = requests.get(url)
        response.raise_for_status()

        outfile = astropy.utils.data.download_file(url, cache=True,
            show_progress=True)

        if outfile!=output_path:
            os.rename(outfile, output_path)

        if os.path.isfile(output_path):
            print(f"Successfully downloaded {url} to {output_path}")
            return output_path

    except Exception as e:
        print(f"Error downloading file at {url}:\n{e}")
        return None
