from neospy.wise import (
    MISSION_PHASES,
    WISE_phase_from_jd,
    WISE_phase_from_scan,
)


def test_WISE_phase_from_jd():
    for phase in MISSION_PHASES.values():
        p0 = WISE_phase_from_jd(phase.jd_start)
        assert p0 == phase, p0.name + "  " + phase.name
        p1 = WISE_phase_from_jd(phase.jd_end - 1)
        assert p1 == phase, p1.name + "  " + phase.name

    assert WISE_phase_from_jd(1000) is None


def test_WISE_phase_from_scan():
    for letter in "rst":
        scan_id = "10000" + letter
        assert WISE_phase_from_scan(scan_id) == MISSION_PHASES["Reactivation_2022"]
        scan_id = "45687" + letter
        assert WISE_phase_from_scan(scan_id) is None

    scan_id = "01000a"
    assert WISE_phase_from_scan(scan_id) == MISSION_PHASES["Cryo"]
    scan_id = "08000k"
    assert WISE_phase_from_scan(scan_id) == MISSION_PHASES["3-Band"]
    scan_id = "09000a"
    assert WISE_phase_from_scan(scan_id) == MISSION_PHASES["Post-Cryo"]
    scan_id = "49000a"
    assert WISE_phase_from_scan(scan_id) == MISSION_PHASES["Reactivation_2014"]
