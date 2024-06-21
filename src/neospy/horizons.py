import requests
import os
import numpy as np

from ._core import HorizonsProperties, Covariance

__all__ = ["HorizonsProperties"]


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

        for neospy_v, jpl_v in lookup.items():
            val = [
                float(p["value"])
                for p in props["orbit"]["elements"]
                if p["label"] == jpl_v
            ]
            if len(val) == 0:
                continue
            phys[neospy_v] = val[0]
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
        for neospy_v, jpl_v in lookup.items():
            val = [p["value"] for p in props["phys_par"] if p["name"] == jpl_v]
            if len(val) == 0:
                continue
            phys[neospy_v] = float(val[0])

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
    from .cache import cache_path
    from .mpc import unpack_designation, pack_designation
    import json

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
def _json(self):
    """
    JSON Representation of the object as queried from horizons.

    This will by default come directly from the cached results, but will query if the
    object is not found.
    """
    return _fetch_json(self.desig, update_cache=False)


HorizonsProperties.fetch = fetch
HorizonsProperties.json = _json
