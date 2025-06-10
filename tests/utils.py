import requests
import re
import os
import astropy.utils.data
import gdown
import bs4
import json
import itertools
import warnings

__all__ = ["download_gdrive_file","make_dummy_wcs","download_github_file"]

class _GoogleDriveFile:
    TYPE_FOLDER = "application/vnd.google-apps.folder"

    def __init__(self, id, name, type):
        self.id = id
        self.name = name
        self.type = type
        self.children = []

    def is_folder(self):
        return self.type == _GoogleDriveFile.TYPE_FOLDER

MAX_NUMBER_FILES = 1000

def download_gdrive_file(
    rel_file_path: str,
    output_dir: str | None = None,
    preserve_gdrive_path: bool = True,
    use_cached : bool = True
) -> str | None:
    """
    Download a file from the POTPyRI google drive folder.
    
    Args:
        rel_file_path (str): Path to the file in the POTPyRI google drive repository.
        output_dir (str, optional): Local directory to save the file (default: current directory).
        preserve_gdrive_path (bool, optional): Preserve the relative path in the output directory (default: True)
        shared_folder_url (str, optional): Base URL of the HISPEC_PARVI_Test_Data repository
            Defaults to: "https://drive.google.com/drive/folders/1bxzsnDmCYYTWIvqdK6si0Uj9GXqQjgdn"
        use_cached (bool, optional): Use cached file if it exists (default: True).
    
    Returns:
        str | None: Path to the downloaded file or None if download failed.
    """
    shared_folder_url = "https://drive.google.com/drive/folders/16gC62idNnVgkVVNHBYUnGRjbNAyabwuT"

    # Output dir
    if output_dir is None:
        output_dir = astropy.utils.data._get_download_cache_loc()
    
    # Set the output path
    if preserve_gdrive_path:
        output_path = os.path.join(output_dir, rel_file_path)
    else:
        output_path = os.path.join(output_dir, os.path.basename(rel_file_path))

    if os.path.exists(output_path) and use_cached:
        return output_path
    
    # Ensure the output directory exists
    if output_dir:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
    try:
        file_id = resolve_relative_path(shared_folder_url, rel_file_path)
        url = f"https://drive.google.com/uc?id={file_id}"
        gdown.download(url, output_path, quiet=False)
    except:
        print(f"Failed to download {rel_file_path} from Google Drive.")
        raise
    
    print(f"Successfully downloaded {rel_file_path} to {output_path}")
    return output_path


def resolve_relative_path(shared_folder_url, relative_path):

    def _parse_google_drive_file(url, content):
        folder_soup = bs4.BeautifulSoup(content, features="html.parser")
        encoded_data = None
        for script in folder_soup.select("script"):
            inner_html = script.decode_contents()
            if "_DRIVE_ivd" in inner_html:
                regex_iter = re.compile(r"'((?:[^'\\]|\\.)*)'").finditer(inner_html)
                try:
                    encoded_data = next(itertools.islice(regex_iter, 1, None)).group(1)
                except StopIteration:
                    raise RuntimeError("Couldn't find the folder encoded JS string")
                break
        if encoded_data is None:
            raise RuntimeError("Cannot retrieve folder info")

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            decoded = encoded_data.encode("utf-8").decode("unicode_escape")
        folder_arr = json.loads(decoded)
        folder_contents = [] if folder_arr[0] is None else folder_arr[0]

        sep = " - "
        splitted = folder_soup.title.contents[0].split(sep)
        if len(splitted) >= 2:
            name = sep.join(splitted[:-1])
        else:
            raise RuntimeError("Cannot extract name")

        gdrive_file = _GoogleDriveFile(
            id=url.split("/")[-1],
            name=name,
            type=_GoogleDriveFile.TYPE_FOLDER,
        )

        id_name_type_iter = [
            (e[0], e[2].encode("raw_unicode_escape").decode("utf-8"), e[3])
            for e in folder_contents
        ]
        return gdrive_file, id_name_type_iter

    def _download_and_parse_google_drive_link(sess, url):
        for _ in range(2):
            if "?" in url:
                url += "&hl=en"
            else:
                url += "?hl=en"
            res = sess.get(url)
            if res.status_code != 200:
                return None
            gdrive_file, id_name_type_iter = _parse_google_drive_file(url, res.text)
            for child_id, child_name, child_type in id_name_type_iter:
                child = _GoogleDriveFile(child_id, child_name, child_type)
                if child.is_folder():
                    sub_url = "https://drive.google.com/drive/folders/" + child_id
                    sub_file = _download_and_parse_google_drive_link(sess, sub_url)
                    if sub_file:
                        child.children = sub_file.children
                gdrive_file.children.append(child)
            return gdrive_file
        return None

    def _get_directory_structure(gdrive_file, previous_path):
        directory_structure = []
        for file in gdrive_file.children:
            file.name = file.name.replace(os.sep, "_")
            if file.is_folder():
                directory_structure.append((None, os.path.join(previous_path, file.name)))
                for i in _get_directory_structure(file, os.path.join(previous_path, file.name)):
                    directory_structure.append(i)
            elif not file.children:
                directory_structure.append((file.id, os.path.join(previous_path, file.name)))
        return directory_structure

    with requests.Session() as sess:
        root = _download_and_parse_google_drive_link(sess, shared_folder_url)
        if not root:
            raise RuntimeError("Could not parse shared folder")

        structure = _get_directory_structure(root, root.name)
        base_prefix = structure[0][1].split('/')[0]
        structure_relative = [(fid, path[len(base_prefix) + 1:]) for fid, path in structure]
        match = [entry for entry in structure_relative if entry[1] == relative_path]
        if not match:
            raise FileNotFoundError(f"{relative_path} not found in folder.")
        return match[0][0]  # return file id

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
