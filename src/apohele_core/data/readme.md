Collection of data files which are used by apohele_core.

Files ending in BSP or BPC are are packed binary SPICE kernels.
These are loaded and packed into the final rust binary, and are included automatically.

`mpc_obs.tsv` - Tab separated table containing known ground observatory codes and
    locations. This is also loaded and compiled into the final binary.

`naif_ids.csv` - CSV containing a table of many of the known NAIF identifiers to full
    names. This is not 100% guaranteed to match, as JPL may decide to renumber objects.
    This is unlikely, but possible.

`Leap_Second.dat` - Leap seconds as provided by IERS, downloaded from:
    https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat
    https://data.iana.org/time-zones/data/leap-seconds.list