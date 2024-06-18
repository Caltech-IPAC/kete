Collection of data files which are used by neospy_core.

Files ending in BSP or BPC are are packed binary SPICE kernels.
These are loaded and packed into the final rust binary, and are included automatically.

`mpc_obs.tsv` - Tab separated table containing known ground observatory codes and
    locations. This is also loaded and compiled into the final binary.

`naif_ids.csv` - CSV containing a table of many of the known NAIF identifiers to full
    names. This is not 100% guaranteed to match, as JPL may decide to renumber objects.
    This is unlikely, but possible.

`simult_state_v0.2.1.bin` - A saved SimultaneousState object from v0.2.1, the last
    breaking change in saving. This is used to ensure backward compatibility when
    loading state files.