BANDS: list[float] = [4700.0, 8000.0]
"""Effective wavelength of the NC1 and NC2 bands in nm."""

FOV_WIDTH: float = 7.10
"""Expected effective field of view width in degrees"""

FOV_HEIGHT: float = 1.68
"""Expected effective field of view height in degrees"""

ZERO_MAGS: list[float] = [170.662, 64.13]
"""Zero point magnitude for nc1 and nc2"""

COLOR_CORR = [
    # Tbb    nc1      nc2
    [200.0, 1.50001, 1.17925],
    [220.0, 1.38320, 1.10277],
    [240.0, 1.29964, 1.05175],
    [260.0, 1.23789, 1.01713],
    [280.0, 1.19106, 0.99348],
    [300.0, 1.15477, 0.97740],
    [400.0, 1.05708, 0.95284],
    [500.0, 1.01897, 0.96149],
]
"""Expected color correction required for black body sources"""
