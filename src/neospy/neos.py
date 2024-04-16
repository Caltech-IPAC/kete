from __future__ import annotations
import numpy as np

# Effective wavelength of the NC1 and NC2 bands in nm.
BANDS: list[float] = [4700.0, 8000.0]

FOV_WIDTH: float = 7.10
FOV_HEIGHT: float = 1.68

# Zero point magnitude for nc1 and nc2
ZERO_MAGS: list[float] = [170.662, 64.13]

COLOR_CORR = np.array(
    [  # Tbb    nc1      nc2
        [200.0, 1.50001, 1.17925],
        [220.0, 1.38320, 1.10277],
        [240.0, 1.29964, 1.05175],
        [260.0, 1.23789, 1.01713],
        [280.0, 1.19106, 0.99348],
        [300.0, 1.15477, 0.97740],
        [400.0, 1.05708, 0.95284],
        [500.0, 1.01897, 0.96149],
    ]
)
