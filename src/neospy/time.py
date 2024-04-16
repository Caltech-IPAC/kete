from __future__ import annotations
import datetime
from astropy.time import Time as AstroTime  # type: ignore
import warnings

__all__ = ["Time", "float_day_to_d_h_m_s", "d_h_m_s_to_float_days"]


class Time:
    """
    A representation of time, always in JD with TDB scaling.

    A wrapper around the `astopy.time.Time` class which enables some convenient
    conversions. Specifically days may be represented by floats, where astopy requires
    ints.

    This also suppresses a number of warnings related to datetime conversion in the
    future. Since datetime doesn't handle leap seconds in the future well, there are a
    number of warnings that get raised in astropy's Time class. However these warnings
    would result in no more than single second errors due to the leap second
    miscalculation. These warnings are suppressed in this class.

    Parameters
    ----------
    in_time:
        Accepts multiple formats of input, including floats, strings, datetime objects,
        see `astopy.time.Time` for more details.
    format:
        Accepts all formats specified by astropy.
    scale:
        Accepts all scales specified by astropy, but all scales are converted to TDB
        immediately.
    """

    def __init__(
        self,
        in_time: (
            str
            | float
            | list[int | float]
            | tuple[int, int, float]
            | tuple[int, int, int]
            | datetime.datetime
        ),
        format: str | None = "jd",
        scale: str = "tdb",
    ):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            if isinstance(in_time, AstroTime):
                self.time = in_time.tdb
            else:
                self.time = AstroTime(in_time, format=format, scale=scale).tdb

    @classmethod
    def from_ymd(cls, year: int, month: int, day: float | int, scale: str = "tdb"):
        """
        Create time object from the Year, Month, and Day.

        Parameters
        ----------
        year:
            The Year, for example `2020`
        month:
            The Month as an integer, 0 = January etc.
        day:
            The day as an integer or float.
        scale:
            Accepts all scales specified by astropy, but all scales are converted to TDB
            immediately.
        """
        # The input format should be (year, month, day)
        day, hour, min, sec = float_day_to_d_h_m_s(day)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            return cls(
                AstroTime(
                    dict(
                        year=int(year),
                        month=int(month),
                        day=day,
                        hour=hour,
                        minute=min,
                        second=sec,
                    ),
                    scale=scale,
                ).tdb
            )

    @classmethod
    def J2000(cls) -> "Time":
        return cls(2451545.0)

    @property
    def jd(self) -> float:
        """Julian Date."""
        return self.time.jd

    @property
    def mjd(self) -> float:
        """Modified Julian Date."""
        return self.time.mjd

    @property
    def utc(self) -> AstroTime:
        """Time in the UTC time format."""
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            out_time = self.time.utc
        return out_time

    @property
    def ymd(self) -> tuple:
        """Return (year, month, day), where day is a float."""
        t = self.time.datetime
        return (
            t.year,
            t.month,
            d_h_m_s_to_float_days(
                d=t.day, h=t.hour, m=t.minute, s=t.second + t.microsecond * 1e-6
            ),
        )

    @property
    def iso(self) -> str:
        """Time in the ISO time format."""
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            return self.time.iso

    def local_time(self, timezone=None) -> str:
        """
        String representation of the Time in a localized time zone format.
        This will automatically take into account daylight savings time if necessary.

        Parameters
        ----------
        timezone:
            Optional, a ``datetime.timezone`` object which defines a time zone, if you
            are using python 3.9 or greater, then the ``zoneinfo`` package in base
            python maintains a list of pre-defined timezones:

                from zoneinfo import available_timezones, ZoneInfo
                available_timezones() # List all available timezones
                timezone = ZoneInfo("US/Pacific")
                neospy.Time.J2000().local_time(timezone)

        """

        if timezone is None:
            timezone = datetime.datetime.now().astimezone().tzinfo
        t = self.utc.datetime
        return (t + timezone.utcoffset(t)).strftime("%Y-%m-%d %X.%f")

    @classmethod
    def from_current_time(cls) -> "Time":
        """Create new Time from the current UTC time."""
        return cls(AstroTime.now(), format=None, scale="utc")

    @classmethod
    def from_strptime(cls, date_string, fmt: str, scale: str = "tdb") -> "Time":
        """
        Create a new AstTime from a formatted text string and the format definition, see
        https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes
        for allowed string formats.
        """
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            return cls(
                AstroTime.strptime(date_string, fmt, scale=scale),
                format=None,
                scale=scale,
            )

    def strftime(self, fmt: str) -> str:
        """
        Generate a string representation of this object using using `datetime` string of
        time, see:
        https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes
        for acceptable input string formats.

        For example: `'%Y-%b-%d %H:%M:%S.f %Z'`
        """
        return self.time.strftime(fmt)

    def __repr__(self):
        return f"{type(self).__name__}({self.jd})"


def float_day_to_d_h_m_s(day: float) -> tuple[int, int, int, float]:
    """
    Convert a float of days into ints of (days, hours, minutes, seconds).
    """
    seconds = datetime.timedelta(days=float(day)).total_seconds()
    days = seconds // (24 * 60 * 60)
    rem_sec = seconds % (24 * 60 * 60)
    hours = rem_sec // (60 * 60)
    rem_sec = rem_sec % (60 * 60)
    minutes = rem_sec // 60
    return (int(days), int(hours), int(minutes), rem_sec % 60)


def d_h_m_s_to_float_days(
    d: float = 0, h: float = 0, m: float = 0, s: float = 0
) -> float:
    """
    Convert floats of days, hours, minutes, and seconds into a single float in units of
    days.
    """
    dt = datetime.timedelta(days=d, hours=h, minutes=m, seconds=s).total_seconds()
    return dt / 60 / 60 / 24
