from __future__ import annotations
import requests
import requests.adapters
import zlib
import json
import logging
import urllib
import os
import zipfile
from importlib import resources
from urllib3.util.retry import Retry


CACHE_DIR = os.getenv("NEOSPY_CACHE_DIR", "~/.neospy/")
"""
Cache directory location is set by this global variable.

By default this puts cache information into the home directory under a hidden
folder.

This may be either set by the environmental variable `NEOSPY_CACHE_DIR` or
by changing this variable after importing neospy.
"""

logger = logging.getLogger(__name__)

__all__ = [
    "data_path",
    "cache_path",
    "cached_file_download",
    "cached_gzip_json_download",
    "cached_zip_download",
    "data_ls",
    "get_cache_size",
]


def cache_path(subfolder=""):
    """
    Helper function to return the absolute location of the cache folder.

    The cache folder contains files which are downloaded during use of neospy and
    are not required for basic function.
    """
    path = os.path.abspath(os.path.expanduser(os.path.expandvars(CACHE_DIR)))
    path = os.path.join(path, subfolder)
    if not os.path.isdir(path):
        logger.info("Cache folder does not exist, creating it now. (%s)", path)
        os.makedirs(path)
    return path


def data_path(subfolder="", module_path="neospy.data"):
    """
    Helper function to the absolute path of the neospy data path.

    The data path is where core data files are stored for neospy, such as required
    spice kernels.
    """
    # Python > 3.9 needs this
    try:
        return os.path.abspath(resources.files(module_path).joinpath(subfolder))
    except AttributeError:
        pass

    # python <= 3.9 needs this
    with resources.path(module_path, "") as p:
        return os.path.abspath(p / subfolder)


def cache_ls(subfolder=""):
    """List the contents of the neospy cache directory"""
    path = os.path.join(cache_path(), subfolder)
    if not os.path.isdir(path):
        return ValueError(f"({path} is not a directory.)")
    return sorted(os.listdir(path))


def data_ls():
    """List the contents of the neospy data directory specified."""
    return sorted(os.listdir(data_path(".")))


def get_cache_size():
    """
    Return the total size of the cache directory in MB.
    """
    return (
        sum(
            os.path.getsize(os.path.join(dirpath, filename))
            for dirpath, dirnames, filenames in os.walk(cache_path())
            for filename in filenames
        )
        / 1024**2
    )


def cached_file_download(url, force_download=False, subfolder=""):
    """
    Download a file from the specified URL.

    This operation is cached, so requesting the same URL will result in the previously
    fetched results being returned. Setting force_download to true will force the cached
    file to be re-downloaded.
    """
    parsed = urllib.parse.urlparse(url)
    filename = os.path.join(cache_path(subfolder), parsed.path.split("/")[-1])
    # check if already downloaded
    if not os.path.isfile(filename) or force_download:
        with requests.Session() as session:
            retry = Retry(connect=3, backoff_factor=0.5)
            adapter = requests.adapters.HTTPAdapter(max_retries=retry)
            session.mount("http://", adapter)
            session.mount("https://", adapter)
            with session.get(url, stream=True, timeout=10) as res:
                res.raise_for_status()
                logger.info("Downloading file from (%s)", url)

                with open(filename, "wb") as f:
                    for chunk in res.iter_content(chunk_size=None):
                        if chunk:
                            f.write(chunk)
    else:
        logger.debug("Previously cached file (%s)", filename)
    return filename


def cached_gzip_json_download(url, force_download=False, subfolder=""):
    """
    Download a gzipped json file from the specified URL.

    This operation is cached, so requesting the same URL will result in the previously
    fetched results being returned. Setting force_download to true will force the cached
    file to be re-downloaded.
    """
    filename = cached_file_download(url, force_download, subfolder)
    # unpack the gzip, then the json
    with open(filename, "rb") as f:
        raw_data = zlib.decompress(f.read(), 16 + zlib.MAX_WBITS).decode()
    return json.loads(raw_data)


def cached_zip_download(url, force_download=False, subfolder=""):
    """
    Download a zipped file to the `neospy/data` directory and unzip it in place.
    Returning the path to the folder where it is unzipped.

    This operation is cached, so requesting the same URL will result in the previously
    fetched results being returned. Setting force_download to true will force the cached
    file to be re-downloaded.
    """
    filename = cached_file_download(url, force_download, subfolder)
    parent_dir, file = os.path.split(filename)
    zipfile.ZipFile(filename).extractall(parent_dir)
    return os.path.join(parent_dir, file.replace(".zip", ""))
