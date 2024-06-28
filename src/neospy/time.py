"""
Time conversion support, primarily provided by the :py:class:`Time` class.
"""

from __future__ import annotations
import calendar
import datetime
import warnings
from astropy.time import Time as AstroTime

__all__ = ["Time", "float_day_to_d_h_m_s", "d_h_m_s_to_float_days", "days_in_year"]


class Time:
    """
    A representation of time, always in JD with TDB scaling.

    A wrapper around the `astropy.time.Time` class which enables some convenient
    conversions. Specifically days may be represented by floats, where Astropy requires
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
        see `astropy.time.Time` for more details.
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
        format: str | None = "jd",  # pylint: disable=redefined-builtin
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
        day, hour, minute, sec = float_day_to_d_h_m_s(day)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            return cls(
                AstroTime(
                    dict(
                        year=int(year),
                        month=int(month),
                        day=day,
                        hour=hour,
                        minute=minute,
                        second=sec,
                    ),
                    scale=scale,
                ).tdb
            )

    @classmethod
    def J2000(cls) -> "Time":
        """
        J2000 is typically defined in TT time scaling, which is ~4.8 milliseconds
        different than TDB.
        """
        return cls(2451545.0, scale="tt")

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
        """
        Return (year, month, day), where day is a float.

        >>> neospy.Time.from_ymd(2010, 1, 1).ymd
        (2010, 1, 1.0)
        """
        t = self.time.datetime
        return (
            t.year,
            t.month,
            d_h_m_s_to_float_days(
                days=t.day,
                hours=t.hour,
                minutes=t.minute,
                seconds=t.second + t.microsecond * 1e-6,
            ),
        )

    @property
    def year_float(self) -> float:
        """
        Time as the UTC year in float form.

        Note that Time is TDB Scaled, causing UTC to be a few seconds different.

        >>> neospy.Time.from_ymd(2010, 1, 1).year_float
        2009.999997875444

        >>> neospy.Time(2457754.5, scale='utc').year_float
        2017.0

        2016 was a leap year, so 366 days instead of 365.

        >>> neospy.Time(2457754.5 - 366, scale='utc').year_float
        2016.0

        """
        tup = self.utc.datetime.timetuple()
        frac = d_h_m_s_to_float_days(
            tup.tm_yday - 1, tup.tm_hour, tup.tm_min, tup.tm_sec
        )
        return tup.tm_year + frac / days_in_year(tup.tm_year)

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

    >>> neospy.time.float_day_to_d_h_m_s(1.5)
    (1, 12, 0, 0.0)

    """
    seconds = datetime.timedelta(days=float(day)).total_seconds()
    days = seconds // (24 * 60 * 60)
    rem_sec = seconds % (24 * 60 * 60)
    hours = rem_sec // (60 * 60)
    rem_sec = rem_sec % (60 * 60)
    minutes = rem_sec // 60
    return (int(days), int(hours), int(minutes), rem_sec % 60)


def d_h_m_s_to_float_days(
    days: float = 0, hours: float = 0, minutes: float = 0, seconds: float = 0
) -> float:
    """
    Convert floats of days, hours, minutes, and seconds into a single float in units of
    days.

    >>> neospy.time.d_h_m_s_to_float_days(1, 11, 59, 60)
    1.5

    """
    dt = datetime.timedelta(
        days=days, hours=hours, minutes=minutes, seconds=seconds
    ).total_seconds()
    return dt / 60 / 60 / 24


def days_in_year(year: int) -> int:
    """
    Return the number of days in the specified year, specifically for leap years.

    >>> neospy.time.days_in_year(2016)
    366

    >>> neospy.time.days_in_year(2017)
    365

    Parameters
    ----------
    year :
        Year as an int.
    """
    return 365 + calendar.isleap(year)
