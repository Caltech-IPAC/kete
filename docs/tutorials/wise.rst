WISE Precovery Calculations
===========================

This example loads the last release year of NEOWISE fields of view, and predicts
which known numbered NEOs are present in all of the frames of NEOWISE during that time.

This can be thought of as an extended version of the KONA tutorial. 


WISE Fields of View
-------------------

A field of view is a patch of sky as seen from an observer. kete supports downloading
these FOVs directly from IPACs IRSA data repository.


.. code-block:: python

    import kete

    # During 2023 there were 2.369 million individual frames taken of both W1 and W2
    # bands, totalling just shy of 5 million individual frames.
    # This may take some time to download from IRSA.
    fovs = kete.wise.fetch_WISE_fovs('Reactivation_2023')

Minor Planet Center Data
------------------------

Now we have to download the orbit data of NEOs from the Minor Planet Center (MPC).


.. code-block:: python

    # This is the orbit table from the MPC, this is over 100Mb, which also may take a
    # minute to download. Note that kete caches these files after the first download.
    # Future queries will use the cached copy unless `force_download=True` is set.
    orb_table = kete.mpc.fetch_known_orbit_data()

    # Now we can down select these to just NEOs:
    is_neo = kete.population.neo(orb_table['peri_dist'], orb_table['ecc'])
    subset = orb_table[is_neo]

    # select only the numbered asteroids
    subset = subset[[str(n).isdigit() for n in subset.desig]]

    # Convert that table of values into cartesian state vectors:
    states = kete.mpc.table_to_states(subset)


Propagation
-----------

Now that we have the states loaded, we have to do a little preparation.

The MPC orbit table will have a mixture of epochs for the objects, we need to
first bring all of the epochs to a common time. A convenient time is the first
frame in the mission phase we have selected.

.. code-block:: python

    # Time of the first exposure.
    jd = fovs[0].jd

    # now we propagate the NEOs to that time, including the effects of the 5 largest
    # main belt asteroids to include more precision. This may take a few minutes.
    states = kete.propagate_n_body(states, jd, include_asteroids=True)

Visibility Test
---------------

Now that we have the FOVs and the States, we can compute what should have been
observable to NEOWISE during this mission phase.

.. note::

    This takes a while to run, as it is computing the State of all NEOs for
    every ~2.4 million frames during 2023 of NEOWISE. On a high end desktop this
    takes around 3 minutes, on a relatively standard macbook pro expect it to
    take around 10 minutes. This time scales relatively linearly with number
    of objects and the number of fields of view.
    
.. code-block:: python

    # This accepts States and a list of field of views, and computes when objects
    # would have been seen. The final output will only include times where objects
    # were seen, and it drops empty FOVs.

    # Compute observable objects.
    visible = kete.fov_state_check(states, fovs)


.. note::

    The outputs of this may be saved using the following:
    
    ``kete.SimultaneousStates.save_list(visible, "visible_wise_2023.bin")``

    The states may later be loaded using:

    ``visible = kete.SimultaneousStates.load_list("visible_wise_2023.bin")``


Computing Positions
-------------------

We can now compute the on-sky positions of these objects as seen from NEOWISE.

Here is a codeblock which prints the first `n_show=100` objects.

.. code-block:: python
        
    n_show = 100
    print("Found: ", len(visible))
    print(f"Displaying the first {n_show}")
    print(f"{'Name':<16}{'mjd':<16}{'RA':<16}{'DEC':<16}{'scan-frame':<16}")
    print("-"*(16 * 5))
    for vis in visible[:n_show]:
        for state in vis:
            vec = (state.pos - vis.fov.observer.pos).as_equatorial
            mjd = kete.Time(vis.fov.jd).mjd
            print((f"{state.desig:<15s},{mjd:<15.6f},{vec.ra_hms:<15s},"
                   f"{vec.dec_dms:<15s},{vis.fov.scan_id}-{str(vis.fov.frame_num)}"))


::

    Found:  68447
    Displaying the first 100
    Name            mjd             RA              DEC             scan-frame      
    --------------------------------------------------------------------------------
    489453         ,59945.005479   ,01 08 21.420   ,+30 49 30.31   ,46370r-175
    279816         ,59945.015411   ,20 22 46.492   ,+69 13 35.87   ,46370r-261
    279816         ,59945.015538   ,20 22 46.540   ,+69 13 35.83   ,46370r-262
    254417         ,59945.016939   ,18 54 08.690   ,+68 51 49.07   ,46370r-274
    162926         ,59945.026871   ,13 45 32.700   ,+31 30 08.87   ,46372r-54
    4544           ,59945.029291   ,13 17 39.888   ,+19 49 46.23   ,46372r-75
    513572         ,59945.030437   ,13 08 11.160   ,+14 19 55.25   ,46372r-85
    455594         ,59945.030819   ,13 05 59.809   ,+12 09 07.12   ,46372r-88
    550271         ,59945.032219   ,12 51 36.202   ,+05 15 39.04   ,46372r-100
    620064         ,59945.032856   ,12 46 41.624   ,+01 59 36.27   ,46372r-106
    277810         ,59945.035403   ,12 25 45.537   ,-10 45 33.28   ,46372r-128
    455687         ,59945.064054   ,02 06 58.093   ,-02 02 57.47   ,46373r-93
    506491         ,59945.065964   ,01 54 09.057   ,+07 51 45.45   ,46373r-110
    163373         ,59945.066983   ,01 46 19.075   ,+13 17 38.21   ,46373r-119
    427621         ,59945.066983   ,01 46 29.691   ,+13 21 45.22   ,46373r-119
    416151         ,59945.068002   ,01 37 04.043   ,+17 55 32.28   ,46373r-127
    434633         ,59945.069403   ,01 25 37.857   ,+25 20 07.03   ,46373r-139
    138852         ,59945.069657   ,01 23 19.125   ,+26 28 07.00   ,46373r-142
    279816         ,59945.080608   ,20 22 59.741   ,+69 13 18.27   ,46373r-236
    162926         ,59945.092068   ,13 45 36.279   ,+31 30 46.71   ,46374r-78
    455594         ,59945.095888   ,13 06 07.171   ,+12 07 36.26   ,46374r-111
    455594         ,59945.096016   ,13 06 07.186   ,+12 07 36.18   ,46374r-113
    495833         ,59945.097926   ,12 50 40.550   ,+02 05 59.83   ,46374r-129
    1627           ,59945.098308   ,12 47 13.930   ,+00 13 16.36   ,46374r-132
    277810         ,59945.100473   ,12 25 50.805   ,-10 46 03.55   ,46374r-151
    378842         ,59945.102128   ,12 12 38.845   ,-19 14 51.80   ,46374r-165
    162082         ,59945.104038   ,11 53 32.643   ,-28 42 27.74   ,46374r-182
    8566           ,59945.121611   ,03 21 18.554   ,-39 41 22.55   ,46375r-52
    481918         ,59945.130143   ,02 04 42.515   ,+03 23 43.26   ,46375r-125
    194268         ,59945.130270   ,02 01 54.844   ,+04 00 00.79   ,46375r-126
    162926         ,59945.157138   ,13 45 39.847   ,+31 31 24.48   ,46376r-24
    4544           ,59945.159558   ,13 17 53.127   ,+19 49 19.41   ,46376r-45
    513572         ,59945.160704   ,13 08 34.064   ,+14 20 06.50   ,46376r-55
    455594         ,59945.161086   ,13 06 14.552   ,+12 06 05.19   ,46376r-59
    550271         ,59945.162486   ,12 51 45.625   ,+05 13 36.48   ,46376r-71
    620064         ,59945.163123   ,12 46 44.949   ,+01 59 31.30   ,46376r-76
    277810         ,59945.165670   ,12 25 56.080   ,-10 46 33.69   ,46376r-98
    481918         ,59945.195340   ,02 04 41.879   ,+03 23 44.14   ,46377r-36
    162926         ,59945.222335   ,13 45 43.417   ,+31 32 02.40   ,46378r-77
    455594         ,59945.226155   ,13 06 21.922   ,+12 04 34.05   ,46378r-110
    455594         ,59945.226283   ,13 06 21.937   ,+12 04 33.97   ,46378r-112
    495833         ,59945.228193   ,12 50 44.136   ,+02 05 47.78   ,46378r-128
    1627           ,59945.228575   ,12 47 23.547   ,+00 12 49.66   ,46378r-131
    277810         ,59945.230740   ,12 26 01.336   ,-10 47 03.85   ,46378r-150
    378842         ,59945.232395   ,12 12 54.299   ,-19 17 15.88   ,46378r-164
    162082         ,59945.234305   ,11 53 44.379   ,-28 44 43.02   ,46378r-181
    8566           ,59945.251878   ,03 21 18.030   ,-39 38 24.29   ,46379r-51
    482650         ,59945.254807   ,02 52 08.786   ,-24 45 03.47   ,46379r-76
    530531         ,59945.258627   ,02 22 16.836   ,-05 39 34.46   ,46379r-109
    497230         ,59945.260664   ,02 06 51.136   ,+04 54 12.69   ,46379r-127
    441641         ,59945.262192   ,01 54 24.948   ,+12 50 35.97   ,46379r-140
    475950         ,59945.262192   ,01 55 21.476   ,+12 40 29.83   ,46379r-140
    424392         ,59945.262829   ,01 48 55.332   ,+16 24 59.36   ,46379r-145
    424392         ,59945.262956   ,01 48 55.339   ,+16 24 59.20   ,46379r-147
    254417         ,59945.277472   ,18 55 02.214   ,+68 58 08.66   ,46379r-272
    162926         ,59945.287405   ,13 45 46.976   ,+31 32 40.25   ,46380r-52
    513572         ,59945.290970   ,13 08 56.906   ,+14 20 17.84   ,46380r-83
    455594         ,59945.291353   ,13 06 29.312   ,+12 03 02.68   ,46380r-86
    550271         ,59945.292753   ,12 51 55.028   ,+05 11 33.88   ,46380r-98
    620064         ,59945.293390   ,12 46 48.260   ,+01 59 26.44   ,46380r-104
    277810         ,59945.295937   ,12 26 06.598   ,-10 47 33.89   ,46380r-126
    8566           ,59945.317075   ,03 21 17.785   ,-39 36 55.06   ,46381r-27
    481918         ,59945.325479   ,02 04 40.619   ,+03 23 46.02   ,46381r-100
    481918         ,59945.325607   ,02 04 40.618   ,+03 23 45.99   ,46381r-101
    162926         ,59945.352475   ,13 45 50.530   ,+31 33 18.15   ,46382r-76
    162926         ,59945.352602   ,13 45 50.537   ,+31 33 18.26   ,46382r-78
    455594         ,59945.356422   ,13 06 36.691   ,+12 01 31.25   ,46382r-111
    455594         ,59945.356550   ,13 06 36.706   ,+12 01 31.17   ,46382r-112
    495833         ,59945.358460   ,12 50 47.706   ,+02 05 35.85   ,46382r-128
    1627           ,59945.358842   ,12 47 33.156   ,+00 12 23.09   ,46382r-131
    378842         ,59945.362662   ,12 13 09.773   ,-19 19 40.09   ,46382r-164
    162082         ,59945.364572   ,11 53 56.094   ,-28 46 58.13   ,46382r-181
    8566           ,59945.382145   ,03 21 17.549   ,-39 35 25.90   ,46383r-51
    482650         ,59945.385074   ,02 52 04.871   ,-24 42 54.37   ,46383r-76
    486607         ,59945.388384   ,02 26 09.410   ,-08 05 50.81   ,46383r-105
    530531         ,59945.388894   ,02 22 16.732   ,-05 33 49.62   ,46383r-109
    497230         ,59945.390931   ,02 06 50.606   ,+04 54 34.00   ,46383r-127
    441641         ,59945.392459   ,01 54 23.878   ,+12 50 51.65   ,46383r-140
    475950         ,59945.392459   ,01 55 21.502   ,+12 40 40.27   ,46383r-140
    424392         ,59945.393096   ,01 49 00.171   ,+16 25 21.89   ,46383r-145
    424392         ,59945.393223   ,01 49 00.177   ,+16 25 21.74   ,46383r-147
    199003         ,59945.399081   ,00 46 38.701   ,+46 06 53.28   ,46383r-197
    254417         ,59945.407613   ,18 55 29.088   ,+69 01 18.50   ,46383r-271
    254417         ,59945.407740   ,18 55 29.130   ,+69 01 18.73   ,46383r-272
    162926         ,59945.417672   ,13 45 54.088   ,+31 33 56.20   ,46384r-52
    513572         ,59945.421110   ,13 09 19.660   ,+14 20 29.01   ,46384r-82
    513572         ,59945.421237   ,13 09 19.683   ,+14 20 29.27   ,46384r-83
    455594         ,59945.421619   ,13 06 44.089   ,+11 59 59.60   ,46384r-86
    277810         ,59945.426204   ,12 26 17.093   ,-10 48 33.90   ,46384r-126
    8566           ,59945.447342   ,03 21 17.325   ,-39 33 56.61   ,46385r-27
    162926         ,59945.482742   ,13 45 57.634   ,+31 34 34.18   ,46386r-77
    455594         ,59945.486689   ,13 06 51.476   ,+11 58 27.88   ,46386r-111
    495833         ,59945.488727   ,12 50 51.260   ,+02 05 24.03   ,46386r-128
    1627           ,59945.489109   ,12 47 42.755   ,+00 11 56.66   ,46386r-132
    378842         ,59945.492929   ,12 13 25.265   ,-19 22 04.43   ,46386r-165
    162082         ,59945.494839   ,11 54 07.788   ,-28 49 13.07   ,46386r-181
    486607         ,59945.518651   ,02 26 08.934   ,-08 04 56.17   ,46387s-16
    530531         ,59945.519161   ,02 22 16.691   ,-05 28 04.66   ,46387s-20
    497230         ,59945.521198   ,02 06 50.096   ,+04 54 55.41   ,46387s-38
    441641         ,59945.522726   ,01 54 22.847   ,+12 51 07.51   ,46387s-51
    475950         ,59945.522726   ,01 55 21.552   ,+12 40 50.83   ,46387s-51
    424392         ,59945.523363   ,01 49 05.060   ,+16 25 44.67   ,46387s-57
    138846         ,59945.525146   ,01 34 52.520   ,+25 36 59.70   ,46387s-72
    199003         ,59945.529348   ,00 47 17.747   ,+45 52 21.25   ,46387s-108