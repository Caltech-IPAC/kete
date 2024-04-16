import pytest
import neospy


@pytest.fixture(scope="session")
def ceres_properties():
    """HorizonsObjectProperties built from Ceres"""
    return neospy.HorizonsProperties(
        desig="Ceres",
        group=None,
        vis_albedo=0.09,
        diameter=939.4,
        moid=0.0,
        peri_dist=2.543224690,
        eccentricity=0.079306796,
        inclination=10.585143108,
        lon_of_ascending=80.15464972,
        peri_arg=72.08042470706,
        peri_time=2463270.778809300,
        h_mag=3.33,
        g_phase=0.12,
        epoch=2462583.0,
        arc_len=None,
        covariance=None,
    )


def test_object_properties(ceres_properties):
    ceres_properties.state
    assert ceres_properties.moid == 0.0


@pytest.mark.horizons
def test_fetch_horizons():
    obj = neospy.HorizonsProperties.fetch("Ceres")
    _ = obj.state
    _ = obj.elements
