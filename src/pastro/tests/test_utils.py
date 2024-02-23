from ..types import Date, CartesianState, KeplerianState
from ..utils import date2jd, jd2date, cartesian2keplerian, keplerian2cartesian
from sgp4.api import jday
import numpy as np


def test_Date():
    """Tests for Date dataclass"""

    # Test tolist method
    date = Date(22, 6, 2019, 16, 12, 56)
    assert date.tolist() == [2019, 6, 22, 16, 12, 56]


def test_date2jd():
    """Test date2jd function."""

    # Test date and expected result
    date = Date(22, 6, 2019, 16, 12, 56)

    # Expected result from sgp4.api.jday
    jd, fr = jday(date.year, date.month, date.day,
                  date.hour, date.minute, date.second)

    # Test function
    assert date2jd(date)[0] == jd + fr
    assert date2jd(date, frac=True) == (jd, fr)
    assert date2jd(date, frac=False)[0] == jd + fr
    assert date2jd(date, mjd=True, frac=False)[0] == jd - 2400000.5 + fr
    assert date2jd(date, mjd=True, frac=True) == (jd - 2400000.5, fr)


def test_jd2date():
    """Test jd2date function."""

    # Test date and expected result
    date = Date(22, 6, 2019, 16, 12, 58)
    jd, fr = jday(date.year, date.month, date.day,
                  date.hour, date.minute, date.second)

    # Test function
    assert jd2date(jd, fr, mjd=False)[0] == date
    assert jd2date(jd - 2400000.5, fr, mjd=True)[0] == date
    assert jd2date(jd + fr, None, mjd=False)[0] == date
    assert jd2date(jd - 2400000.5 + fr, None, mjd=True)[0] == date

    # Test function with error estimate
    _date, eps = jd2date(jd, fr, mjd=False)
    assert _date == date
    assert eps < 1.


def test_keplerian2cartesian():

    mu_sun = 1.3271276489987138E+20

    elements = KeplerianState(1.082073518745249E+11, 6.737204890591715E-03,
                              3.394387859410573E+00, 7.661212333458995E+01,
                              5.520309419058912E+01, 3.559064425112629E+02)

    state = CartesianState(-6.561639868572587E+10, 8.498200549242477E+10,
                           4.953209922912188E+09, -2.784002901222631E+04,
                           -2.159352048214584E+04, 1.309840720051276E+03)

    _state = keplerian2cartesian(elements, mu_sun)
    _elements = cartesian2keplerian(state, mu_sun)

    # Relative errors
    TOL = 2e-13
    assert abs(_state.x - state.x) / np.abs(state.x) < TOL
    assert abs(_state.y - state.y) / np.abs(state.y) < TOL
    assert abs(_state.z - state.z) / np.abs(state.z) < TOL
    assert abs(_state.dx - state.dx) / np.abs(state.dx) < TOL
    assert abs(_state.dy - state.dy) / np.abs(state.dy) < TOL
    assert abs(_state.dz - state.dz) / np.abs(state.dz) < TOL

    assert abs(_elements.a - elements.a) / elements.a < TOL
    assert abs(_elements.e - elements.e) / elements.e < TOL
    assert abs(_elements.i - elements.i) / elements.i < TOL
    assert abs(_elements.Omega - elements.Omega) / elements.Omega < TOL
    assert abs(_elements.omega - elements.omega) / elements.omega < TOL
    assert abs(_elements.nu - elements.nu) / elements.nu < TOL
