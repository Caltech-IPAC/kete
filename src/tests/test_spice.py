import numpy as np
import pytest

from neospy import SpiceKernels, State, Time

# Turn off flake8 checker for this file.
# flake8: noqa

# fmt: off
# Here are expected interpolation equatorial state values as calculated by Horizon.
# [target, jd, pos.x, pos.y, pos.z,
#              vel.x, vel.y, vel.z]
EXPECTED_EQUITORIAL = [
    ["Earth" , 2505035.5,   -0.143102566296,    -0.922493350579,    -0.399580836275,
                             0.016758306839,    -0.002280675804,    -0.000989096637],
    ["Moon"  , 2505035.5,-0.145246016114196,-0.9209739949551414,-0.3990390962430936,
                          0.016436330755163,-0.0027213958479751,-0.0011310662345952],
    ["earth ba", 2485075.5, 0.890691383427617, 0.4087920914165625, 0.177113518664760,
                         -0.007976016046565, 0.0140578719519273, 0.0060915119230729],
    ["Mars"  , 2435040.5, 1.357638349280660,-0.2305928688109060,-0.1426082027542960,
                          0.003255868695805, 0.0135970078786732, 0.0061478371780576],
]
# fmt: on


class TestSpiceKernels:
    @pytest.mark.parametrize("expected", EXPECTED_EQUITORIAL)
    def test_ecliptic_state(self, expected):
        e_pos = expected[2:5]
        e_vel = expected[5:8]
        state = SpiceKernels.state(expected[0], expected[1])
        assert isinstance(state, State)
        assert np.allclose(state.pos.as_equatorial, e_pos)
        assert np.allclose(state.vel.as_equatorial, e_vel)

        # Ensure that the same values are calculated if using an astropy time
        jd = Time(expected[1], format="jd", scale="tdb")
        state = SpiceKernels.state(expected[0], jd)
        assert isinstance(state, State)
        assert np.allclose(state.pos.as_equatorial, e_pos)
        assert np.allclose(state.vel.as_equatorial, e_vel)

    def test_loaded_objects(self):
        assert len(SpiceKernels.loaded_objects()) > 8

    def test_earth_to_ecliptic(self):
        state = SpiceKernels.earth_pos_to_ecliptic(2454832.5, 0, 0, 0, center="399")
        assert np.allclose(
            list(-state.pos),
            [7.686364937857906e-06, -3.847842319488030e-05, 1.667609317250255e-05],
        )
        assert np.allclose(
            -state.vel,
            [2.642169198371849e-04, 4.433586059319487e-05, -1.948260054980404e-05],
        )
