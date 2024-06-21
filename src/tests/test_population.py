import numpy as np
from neospy import population, power_law


def test_population():
    assert population.which_group([1000], [0.0]) == ['distant']
    assert not population.jup_trojan(1, 0)


def test_power_law():
    law = power_law.CumulativePowerLaw([1, 1], [1], 100, 10000)
    assert law.num_larger_than([1]) == [100]
    assert np.allclose(law.value_of_nth_largest([100]),[1])