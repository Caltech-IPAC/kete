import requests
import os
from typing import Union, Optional

import base64
import json
import numpy as np

from ._core import HorizonsProperties, Covariance, NonGravModel, CometElements
from .mpc import unpack_designation, pack_designation
from .time import Time
from .cache import cache_path
from .covariance import generate_sample_from_cov

__all__ = ["HorizonsProperties", "fetch_spice_kernel"]


_PARAM_MAP = {
    "a1": "a1",
    "a2": "a2",
    "a3": "a3",
    "aln": "alpha",
    "nm": "m",
    "r0": "r_0",
    "nk": "k",
    "nn": "n",
    "dt": "dt",
    "e": "eccentricity",
    "q": "peri_dist",
    "tp": "peri_time",
    "node": "lon_of_ascending",
    "peri": "peri_arg",
    "i": "inclination",
}


def fetch(name, update_name=True, cache=True, update_cache=False, exact_name=False):
    """
    Fetch object properties from JPL Horizons service.

    This does a best effort to fill in the values from horizons, but some values
    may be missing. The keyword arguments in this function allow to fill in
    defaults when the specified values cannot be found from horizons.


    Parameters
    ----------
    name :
        Name of the object to fetch.
    update_name :
        If this is true, replace name with the name that horizons uses for the
        object.
    cache :
        If this is true, the results are cached, and no query to horizons takes
        place if the object is already in the local cache.
    update_cache :
        If this is true, the current cache contents are ignored, and new values
        are saved after querying horizons.
    exact_name :
        If this is true, it is assumed that an exact designation in the format of
        horizons has been provided. This is particularly useful for comets which
        have fragments, as these objects are difficult to query since there are
        several names which have matches on substrings.
    """
    props = _fetch_json(
        name,
        update_cache=update_cache,
        update_name=update_name,
        cache=cache,
        exact_name=exact_name,
    )
    if update_name:
        horizons_name = props.get("object", {})
        if "des" in horizons_name:
            name = horizons_name["des"]
        elif "fullname" in horizons_name:
            name = horizons_name["fullname"].split("(")[-1].replace(")", "")
    phys = {
        "desig": name,
        "inclination": None,
        "eccentricity": None,
        "lon_of_ascending": None,
        "h_mag": None,
        "g_phase": None,
        "diameter": None,
        "vis_albedo": None,
        "moid": None,
        "peri_dist": None,
        "peri_arg": None,
        "peri_time": None,
        "arc_len": None,
        "epoch": None,
        "covariance": None,
    }
    if "orbit" in props:
        lookup = {
            "eccentricity": "e",
            "inclination": "i",
            "lon_of_ascending": "node",
            "peri_arg": "peri",
            "peri_dist": "q",
            "peri_time": "tp",
        }
        lookup_rev = dict(x[::-1] for x in lookup.items())

        for kete_v, jpl_v in lookup.items():
            val = [
                float(p["value"])
                for p in props["orbit"]["elements"]
                if p["label"] == jpl_v
            ]
            if len(val) == 0:
                continue
            phys[kete_v] = val[0]

        phys["epoch"] = float(props["orbit"]["epoch"])
        if "moid" in props["orbit"] and props["orbit"]["moid"] is not None:
            phys["moid"] = float(props["orbit"]["moid"])
        if "data_arc" in props["orbit"] and props["orbit"]["data_arc"] is not None:
            phys["arc_len"] = float(props["orbit"]["data_arc"])
        if props["orbit"]["covariance"] is not None:
            cov_epoch = float(props["orbit"]["covariance"]["epoch"])
            mat = np.nan_to_num(
                np.array(props["orbit"]["covariance"]["data"], dtype=float)
            )
            labels = props["orbit"]["covariance"]["labels"]
            if "elements" in props["orbit"]["covariance"]:
                elements = {
                    x["label"]: float(x["value"])
                    for x in props["orbit"]["covariance"]["elements"]
                }
            else:
                elements = {
                    lab: phys[lookup_rev[lab]] for lab in labels if lab in lookup_rev
                }
            if "model_pars" in props["orbit"]:
                for param in props["orbit"]["model_pars"]:
                    elements[param["name"]] = float(param["value"])
            params = [(lookup_rev.get(x, x), elements.get(x, np.nan)) for x in labels]
            phys["covariance"] = Covariance(name, cov_epoch, params, mat)
    else:
        raise ValueError(
            "Horizons did not return orbit information for this object:" f"\n{props}"
        )

    if "phys_par" in props:
        lookup = {
            "h_mag": "H",
            "g_phase": "G",
            "diameter": "diameter",
            "vis_albedo": "albedo",
        }
        for kete_v, jpl_v in lookup.items():
            val = [p["value"] for p in props["phys_par"] if p["name"] == jpl_v]
            if len(val) == 0:
                continue
            phys[kete_v] = float(val[0])

    phys["group"] = props.get("object", {}).get("orbit_class", {}).get("name", None)
    return HorizonsProperties(**phys)


def _fetch_json(
    name, update_name=True, cache=True, update_cache=False, exact_name=False
):
    """
    Fetch object properties from JPL Horizons service.

    This does a best effort to fill in the values from horizons, but some values
    may be missing. The keyword arguments in this function allow to fill in
    defaults when the specified values cannot be found from horizons.


    Parameters
    ----------
    name :
        Name of the object to fetch.
    update_name :
        If this is true, replace name with the name that horizons uses for the
        object.
    cache :
        If this is true, the results are cached, and no query to horizons takes
        place if the object is already in the local cache.
    update_cache :
        If this is true, the current cache contents are ignored, and new values
        are saved after querying horizons.
    exact_name :
        If this is true, it is assumed that an exact designation in the format of
        horizons has been provided. This is particularly useful for comets which
        have fragments, as these objects are difficult to query since there are
        several names which have matches on substrings.
    """

    dir_path = os.path.join(cache_path(), "horizons_props")
    try:
        filename = os.path.join(dir_path, f"{pack_designation(name)}.json")
    except (SyntaxError, ValueError):
        filename = os.path.join(dir_path, f"{name.replace('/', '_')}.json")

    if not os.path.isdir(dir_path) and cache:
        os.makedirs(dir_path)

    if os.path.isfile(filename) and cache and not update_cache:
        try:
            with open(filename, "r", encoding="utf8") as f:
                return json.load(f)
        except Exception:  # pylint: disable=broad-except
            pass

    if isinstance(name, str):
        try:
            name = unpack_designation(name)
        except (SyntaxError, ValueError):
            pass
    else:
        name = str(name)
    query = "des" if exact_name else "sstr"
    response = requests.get(
        f"https://ssd-api.jpl.nasa.gov/sbdb.api?{query}={name}"
        "&phys-par=true&full-prec=true&cov=mat",
        timeout=30,
    )

    if response.status_code == 300:
        names = [
            des["pdes"] for des in response.json()["list"] if "-" not in des["pdes"]
        ]
        if len(names) == 1:
            return _fetch_json(
                names[0], update_name, cache, update_cache, exact_name=True
            )

    response.raise_for_status()

    response_json = response.json()

    if cache:
        with open(filename, "w", encoding="utf8") as f:
            json.dump(response_json, f)
    return response_json


@property  # type: ignore
def _json(self) -> dict:
    """
    JSON Representation of the object as queried from horizons.

    This will by default come directly from the cached results, but will query if the
    object is not found.
    """
    return _fetch_json(self.desig, update_cache=False)


def _nongrav_params(self) -> dict:
    # default parameters used by jpl horizons for their non-grav models.
    params = {
        "a1": 0.0,
        "a2": 0.0,
        "a3": 0.0,
        "alpha": 0.1112620426,
        "r_0": 2.808,
        "m": 2.15,
        "n": 5.093,
        "k": 4.6142,
        "dt": 0.0,
    }

    orbit = self.json["orbit"]
    if "model_pars" not in orbit:
        return params

    for vals in orbit["model_pars"]:
        if vals["name"].lower() not in params:
            raise ValueError("Unknown non-grav values: ", vals)
        params[_PARAM_MAP[vals["name"].lower()]] = float(vals["value"])
    return params


@property  # type: ignore
def _nongrav(self) -> NonGravModel:
    """
    The non-gravitational forces model from the values returned from horizons.
    """
    params = _nongrav_params(self)
    return NonGravModel.new_comet(**params)


def _sample(self, n_samples):
    """
    Sample the covariance matrix for this object, returning `n_samples` of states and
    non-gravitational models, which may be used to propagate the orbit in time.

    This uses the full covariance matrix in order to correctly sample the object's orbit
    with all correlations, including correlations with the non-gravitational forces.

    Parameters
    ----------
    n_samples :
        The number of samples to take of the covariance.
    """
    matrix = self.covariance.cov_matrix
    epoch = self.covariance.epoch
    samples = generate_sample_from_cov(n_samples, matrix)

    elem_keywords = [
        "eccentricity",
        "inclination",
        "lon_of_ascending",
        "peri_arg",
        "peri_dist",
        "peri_time",
    ]

    labels = self.json["orbit"]["covariance"]["labels"]
    mapped_label = [_PARAM_MAP[k.lower()] for k in labels]

    best_params = _nongrav_params(self)
    for prop in elem_keywords:
        best_params[prop] = getattr(self, prop)
    states = []
    non_gravs = []
    for sample in samples:
        cur_params = {
            label: best_params[label] + d for label, d in zip(mapped_label, sample)
        }
        elem_params = {x: cur_params.pop(x) for x in elem_keywords}
        state = CometElements(self.desig, epoch, **elem_params).state
        non_grav = NonGravModel.new_comet(**cur_params)
        states.append(state)
        non_gravs.append(non_grav)

    return states, non_gravs


HorizonsProperties.fetch = fetch
HorizonsProperties.json = _json
HorizonsProperties.non_grav = _nongrav
HorizonsProperties.sample = _sample


def fetch_spice_kernel(
    name,
    jd_start: Union[Time, float],
    jd_end: Union[Time, float],
    exact_name: bool = False,
    update_cache: bool = False,
    apparition_year: Optional[int] = None,
):
    """
    Download a SPICE kernel from JPL Horizons and save it directly into the Cache.

    .. code-block:: python

        from kete import horizons, Time
        jd_start = Time.from_ymd(1900, 1, 1)
        jd_end = Time.from_ymd(2100, 1, 1)
        horizons.fetch_spice_kernel("10p", jd_start, jd_end)

    Parameters
    ----------
    name :
        Name or integer id value of the object.
    jd_start:
        Start date of the SPICE kernel to download.
    jd_end:
        End date of the SPICE kernel to download.
    exact_name:
        If the specified name is the exact name in Horizons, this can help for
        comet fragments.
    update_cache:
        If the current state of the cache should be ignored and the file
        re-downloaded.
    apparition_year:
        If the object is a comet, retrieve the orbit fit which is previous to this
        specified year. If this is not provided, then default to the most recent
        epoch of orbit fit. Ex: `apparition_year=1980` will return the closest
        epoch before 1980.
    """

    if not isinstance(jd_start, Time):
        jd_start = Time(jd_start)
    if not isinstance(jd_end, Time):
        jd_end = Time(jd_end)

    if isinstance(name, str):
        try:
            name = unpack_designation(name)
        except (SyntaxError, ValueError):
            pass
    else:
        name = str(name)

    query = "des" if exact_name else "sstr"
    # Name resolution using the sbdb database
    name_dat = requests.get(
        f"https://ssd-api.jpl.nasa.gov/sbdb.api?{query}={name}",
        timeout=30,
    )
    if "object" not in name_dat.json():
        raise ValueError("Failed to find object: ", str(name_dat.json()))
    comet = "c" in name_dat.json()["object"]["kind"].lower()

    if comet and apparition_year is None:
        apparition_year = jd_end.ymd[0]

    spk_id = int(name_dat.json()["object"]["spkid"])

    dir_path = os.path.join(cache_path(), "kernels")

    if apparition_year is not None:
        filename = os.path.join(dir_path, f"{spk_id}_epoch_{apparition_year}.bsp")
    else:
        filename = os.path.join(dir_path, f"{spk_id}.bsp")

    if os.path.isfile(filename) and not update_cache:
        return

    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)

    jd_s_str = jd_start.to_datetime().strftime("%Y-%m-%d")
    jd_e_str = jd_end.to_datetime().strftime("%Y-%m-%d")
    cap = f"CAP<{apparition_year}%3B" if comet else ""
    response = requests.get(
        f"https://ssd.jpl.nasa.gov/api/horizons.api?COMMAND='DES={spk_id}%3B{cap}'"
        f"&EPHEM_TYPE=SPK&START_TIME='{jd_s_str}'&STOP_TIME='{jd_e_str}'&CENTER=0",
        timeout=30,
    )

    if response.status_code == 300:
        names = [
            des["pdes"] for des in response.json()["list"] if "-" not in des["pdes"]
        ]
        if len(names) == 1:
            fetch_spice_kernel(names[0], jd_start, jd_end, exact_name=True)

    if response.status_code != 200:
        raise OSError(f"Error from Horizons: {response.json()}")

    if "spk" not in response.json():
        raise ValueError("Failed to fetch file\n:", response.json())

    with open(filename, "wb") as f:
        f.write(base64.b64decode(response.json()["spk"]))
