from neospy.time import Time, d_h_m_s_to_float_days, float_day_to_d_h_m_s


class TestTime:
    def test_init(self):
        t = Time(2460676.5)
        assert t.jd == 2460676.5
        assert t.mjd == 60676.0
        assert t.ymd == (2025, 1, 1)
        assert t.iso == "2025-01-01 00:00:00.000"
        assert t.utc.jd == 2460676.4991992605
        assert Time.from_ymd(2025, 1, 1).jd == t.jd

        assert Time.J2000().jd == 2451545.0
        assert Time.from_current_time().jd > Time.J2000().jd

        s = "%Y-%b-%d %H:%M:%S.%f"
        assert Time.from_strptime(Time.J2000().strftime(s), s).jd == Time.J2000().jd

    def test_d_h_m_s_to_float(self):
        t = d_h_m_s_to_float_days(d=1, h=2, m=3, s=4.5)
        assert float_day_to_d_h_m_s(t) == (1, 2, 3, 4.5)
