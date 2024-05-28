from __future__ import annotations
from typing import Union, Optional
from collections import namedtuple
import base64
import glob
import os
import requests
import numpy as np

from .time import Time
from . import _core  # pylint: disable=no-name-in-module
from .constants import AU_KM
from .cache import cached_file_download, cache_path
from .vector import Frames, State

__all__ = ["SpiceKernels"]

SpkInfo = namedtuple("SpkInfo", "name, jd_start, jd_end, center, frame, spk_type")


def _validate_time(time: Union[float, Time]) -> float:
    """
    Verifies that the time provided is either a `float` or
    :class:`~neospy.time.Time` object.

    Parameters
    ----------
    jd:
        Julian time (TDB) of the desired record.

    Returns
    -------
    float
        The Julian time as a float.

    Raises
    ------
    TypeError
        If the input time is not a `float` or `Time`.
    """

    if isinstance(time, Time):
        return time.jd
    elif isinstance(time, float):
        return time
    try:
        return float(time)
    except Exception as exc:
        raise TypeError("Invalid jd type, use Time or float") from exc


class SpiceKernels:
    """
    Class for estimating the ephemeris for given bodies in the solar system using SPICE
    and the DE440 dataset.

    This allows for the loading of additional spice kernels along with querying.
    """

    _name_cache: dict = {}

    @classmethod
    def state(
        cls,
        target: Union[str, int],
        jd: Union[float, Time],
        center: str = "Sun",
        frame: Frames = Frames.Ecliptic,
    ) -> State:
        """
        Calculates the :class:`~neospy.State` of the target object at the
        specified time `jd`.

        This defaults to the ecliptic heliocentric state, though other centers may be
        chosen.

        Parameters
        ----------
        target:
            The names of the target object, this can include any object name listed in
            :meth:`~neospy.spice.SpiceKernels.loaded_objects`
        jd:
            Julian time (TDB) of the desired record.
        center:
            The center point, this defaults to being heliocentric.
        frame:
            Coordinate frame of the state, defaults to ecliptic.

        Returns
        -------
        State
            Returns the ecliptic state of the target in AU and AU/days.

        Raises
        ------
        ValueError
            If the desired time is outside of the range of the source binary file.
        """
        target, ids = cls.name_lookup(target)
        center, center_id = cls.name_lookup(center)
        jd = _validate_time(jd)
        return _core.spk_state(ids, jd, center_id, frame)

    @classmethod
    def name_lookup(cls, name: Union[int, str]) -> tuple[str, int]:
        """
        Given the provided partial name or integer, find the full name contained within
        the loaded SPICE kernels.

        >>> neospy.SpiceKernels.name_lookup("jupi")
        ('jupiter barycenter', 5)

        >>> neospy.SpiceKernels.name_lookup(10)
        ('sun', 10)

        If there are multiple names, but an exact match, the exact match is returned. In
        the case of ``Earth``, there is also ``Earth Barycenter``, but asking for Earth
        will return the exact match. Putting ``eart`` will raise an exception as there
        are 2 partial matches.

        >>> neospy.SpiceKernels.name_lookup("Earth")
        ('earth', 399)

        >>> neospy.SpiceKernels.name_lookup("Earth b")
        ('earth barycenter', 3)

        Parameters
        ----------
        name :
            Name, partial name, or integer id value of the object.

        Returns
        -------
        tuple :
            Two elements in the tuple, the full name and the integer id value.
        """
        if isinstance(name, str):
            name = name.lower()
        if name in cls._name_cache:
            return cls._name_cache[name]

        # barycenter of the solar system is special
        if name == 0:
            return ("SSB", 0)

        try:
            lookup_name = _core.spk_get_name_from_id(int(name))
        except ValueError:
            lookup_name = name
        lookup_name = lookup_name.lower()

        found = []
        for loaded in cls.loaded_objects():
            loaded_lower = loaded[0].lower()
            # If it is an exact match, finish early
            if lookup_name == loaded_lower:
                cls._name_cache[name] = loaded
                return loaded
            if lookup_name in loaded_lower:
                found.append(loaded)
        found = list(set(found))

        if len(found) == 1:
            cls._name_cache[name] = found[0]
            return found[0]
        elif len(found) > 1:
            raise ValueError(f"Multiple objects match this name {found}")
        raise ValueError(f"No loaded objects which match this name ({name})")

    @staticmethod
    def loaded_objects() -> list[tuple[str, int]]:
        """
        Return the name of all objects which are currently loaded in the SPICE kernels.
        """
        objects = _core.spk_loaded()
        return [(_core.spk_get_name_from_id(o), o) for o in objects]

    @classmethod
    def loaded_object_info(cls, desig: Union[int, str]) -> list[SpkInfo]:
        """
        Return the available SPK information for the target object.

        Parameters
        ----------
        name :
            Name or integer id value of the object.
        """
        name, naif = cls.name_lookup(desig)
        return [SpkInfo(name, *k) for k in _core.spk_available_info(naif)]

    @staticmethod
    def cached_kernel_url_download(url, force_download: bool = False):
        """
        Download the target url into the cache folder of spice kernels.
        """
        cached_file_download(url, force_download=force_download, subfolder="kernels")

    @staticmethod
    def cached_kernel_horizons_download(
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

            from neospy import SpiceKernels, Time
            jd_start = Time.from_ymd(1900, 1, 1)
            jd_end = Time.from_ymd(2100, 1, 1)
            SpiceKernels.cached_kernel_horizons_download("10p", jd_start, jd_end)

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
        from .mpc import unpack_designation

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

        jd_s_str = jd_start.strftime("%Y-%m-%d")
        jd_e_str = jd_end.strftime("%Y-%m-%d")
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
                SpiceKernels.cached_kernel_horizons_download(
                    names[0], jd_start, jd_end, exact_name=True
                )

        if response.status_code != 200:
            raise OSError(f"Error from Horizons: {response.json()}")

        if "spk" not in response.json():
            raise ValueError("Failed to fetch file\n:", response.json())

        with open(filename, "wb") as f:
            f.write(base64.b64decode(response.json()["spk"]))

    @staticmethod
    def cached_kernel_ls():
        """
        List all files contained within the kernels cache folder.
        """
        path = os.path.join(cache_path(), "kernels", "**")
        return glob.glob(path)

    @staticmethod
    def kernel_reload(filenames: list[str], include_cache=False):
        """
        Load the specified spice kernels into memory, this resets the currently loaded
        kernels.

        If `include_cache` is true, this will reload the kernels contained within the
        neospy cache folder as well.

        Parameters
        ----------
        filenames :
            Paths to the specified files to load, this must be a list of filenames.
        include_cache:
            This decides if all of the files contained within the neospy cache should
            be loaded in addition to the specified files.
        """
        _core.spk_reset()
        _core.pck_reset()

        if include_cache:
            cache_files = glob.glob(os.path.join(cache_path(), "kernels", "**.bsp"))
            _core.spk_load(cache_files)

            cache_files = glob.glob(os.path.join(cache_path(), "kernels", "**.bpc"))
            _core.pck_load(cache_files)

        _core.spk_load(filenames)

    @staticmethod
    def spk_header_info(filename: str):
        """
        Return the comments contained within the header of the provided DAF file, this
        includes SPK and PCK files.

        This does not load the contents of the file into memory, it only prints the
        header contents.

        Parameters
        ----------
        filename :
            Path to a DAF file.
        """
        return _core.daf_header_info(filename).replace("\x04", "\n").strip()

    @classmethod
    def mpc_code_to_ecliptic(
        cls, obs_code: str, jd: Union[float, Time], center: str = "Sun", full_name=False
    ) -> State:
        """
        Load an MPC Observatory code as an ecliptic state.

        This only works for ground based observatories.

        Parameters
        ----------
        obs_code:
            MPC observatory code or name of observatory.
        jd:
            Julian time (TDB) of the desired state.
        center:
            The new center point, this defaults to being heliocentric.
        full_name:
            Should the final state include the full name of the observatory or just its
            code.

        Returns
        -------
        State
            Returns the equatorial state of the observatory in AU and AU/days.
        """
        from .mpc import find_obs_code

        obs = find_obs_code(obs_code)
        return cls.earth_pos_to_ecliptic(
            jd,
            geodetic_lat=obs[0],
            geodetic_lon=obs[1],
            height_above_surface=obs[2],
            name=obs[3] if full_name else obs[4],
            center=center,
        )

    @classmethod
    def earth_pos_to_ecliptic(
        cls,
        jd: Union[float, Time],
        geodetic_lat: float,
        geodetic_lon: float,
        height_above_surface: float,
        name: Optional[None] = None,
        center: str = "Sun",
    ) -> State:
        """
        Given a position in the frame of the Earth at a specific time, convert that to
        Sun centered ecliptic state.

        This uses the WGS84 model of Earth's shape to compute state. This uses Geodetic
        latitude and longitude, not geocentric.

        The frame conversion is done using a high precision PCK file from the NAIF/JPL
        website. The file provided in neospy is valid from ~2000 to ~2024. This file
        gives positional error on the scale of a few cm. Additionally a lower resolution
        file is provided which is valid until 2099, and will be used automatically when
        querying past the end of the high resolution file.

        Parameters
        ----------
        jd:
            Julian time (TDB) of the desired state.
        geodetic_lat:
            Latitude on Earth's surface in degrees.
        geodetic_lon:
            Latitude on Earth's surface in degrees.
        height_above_surface:
            Height of the observer above the surface of the Earth in km.
        name :
            Optional name of the position on Earth.
        center:
            The new center point, this defaults to being heliocentric.

        Returns
        -------
        State
            Returns the equatorial state of the target in AU and AU/days.
        """
        jd = _validate_time(jd)
        pos = _core.wgs_lat_lon_to_ecef(
            geodetic_lat, geodetic_lon, height_above_surface
        )
        pos = np.array(pos) / AU_KM
        _, center_id = cls.name_lookup(center)
        return _core.pck_earth_frame_to_ecliptic(pos, jd, center_id, name)

    @classmethod
    def moon_illumination_frac(cls, jd: Union[float, Time]):
        """
        Compute the fraction of the moon which is illuminated at the specified time.

        This is a simple approximation using basic spherical geometry, and assumes
        the observer is located at the geocenter of the Earth.

        >>> SpiceKernels.moon_illumination_frac(Time.from_ymd(2024, 2, 24))
        0.9622388594340928

        Parameters
        ----------
        jd:
            Julian time (TDB) of the desired state.

        Returns
        -------
        State
            Fraction between 0 and 1 of the moons visible surface which is illuminated.
        """
        jd = _validate_time(jd)

        moon2sun = -cls.state("moon", jd).pos
        moon2earth = -cls.state("moon", jd, center="399").pos
        return 1.0 - moon2sun.angle_between(moon2earth) / 180
