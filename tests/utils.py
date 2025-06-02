import requests
import re
import os
from astropy.time import Time
from tqdm import tqdm
import astropy.utils.data

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
    
    # LFS URL
    batch_url = "https://github.com/CIERA-Transients/POTPyRI_test.git/info/lfs/objects/batch"

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

        # Check if the file is an LFS pointer (starts with 'version')
        if response.text.startswith('version'):
            # Extract the oid (SHA-256 hash) and size from the LFS pointer file
            match = re.search(r'oid sha256:([a-f0-9]{64})', response.text)
            size_match = re.search(r'size (\d+)', response.text)
            if not match or not size_match:
                raise ValueError(f"Unable to find 'oid' or 'size' in the LFS pointer file: {response.text}")

            lfs_oid = match.group(1)
            lfs_size = float(size_match.group(1))

            # Build the JSON payload for LFS batch API
            payload = {
                "operation": "download",
                "transfer": ["basic"],
                "objects": [{"oid": lfs_oid, "size": int(lfs_size)}]
            }

            # API URL to interact with LFS objects batch download
            headers = {
                'Accept': 'application/vnd.git-lfs+json',
                'Content-Type': 'application/json',
            }

            # Send the request to GitHub LFS API
            lfs_response = requests.post(lfs_batch_url, headers=headers, json=payload)
            lfs_response.raise_for_status()

            # Extract the URL for the LFS object
            lfs_data = lfs_response.json()
            download_url = lfs_data['objects'][0]['actions']['download']['href']

            # Download the file from the href URL
            #file_data = requests.get(download_url).content
            # Stream the file and use tqdm for progress
            with requests.get(download_url, stream=True) as r, \
                 open(output_path, "wb") as f, \
                 tqdm(total=lfs_size, unit="B", unit_scale=True, desc=f"Downloading {git_file_path}") as pbar:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
                    pbar.update(len(chunk))

        else:

            file_data = response.content
            with open(output_path, 'wb') as f:
               f.write(file_data)
        
        print(f"Successfully downloaded {url} to {output_path}")
        return output_path

    except Exception as e:
        print(f"Error downloading file at {url}:\n{e}")
        return None
