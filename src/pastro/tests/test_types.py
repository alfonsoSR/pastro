from ..types import CartesianState, KeplerianState, SphericalGeocentric
import numpy as np


def test_CartesianState_operations() -> None:

    a = CartesianState(1, 2, 3, 4, 5, 6)
    b = CartesianState(1, 2, 3, 4, 5, 6)

    # Test addition
    assert a + b == CartesianState(2, 4, 6, 8, 10, 12)

    # Test subtraction
    assert a - b == CartesianState(0, 0, 0, 0, 0, 0)

    return None


def test_iterators() -> None:

    x = np.ones((6,))
    cscalar = CartesianState(1, 2, 3, 4, 5, 6)
    cvector = CartesianState(x, 2 * x, 3 * x, 4 * x, 5 * x, 6 * x)

    scalar_count = 0
    for state in cscalar:
        assert isinstance(state, CartesianState)
        assert state == cscalar
        scalar_count += 1

    vector_count = 0
    for state in cvector:
        assert isinstance(state, CartesianState)
        assert state == cscalar
        vector_count += 1

    assert scalar_count == 1
    assert vector_count == 6


def test_CartesianStatesArray_operations() -> None:

    a = CartesianState(np.array([1, 2, 3]), np.array([4, 5, 6]),
                       np.array([0., 1., 7.]), np.array([0., 1., 7.]),
                       np.array([0., 1., 7.]), np.array([0., 1., 7.]))

    b = CartesianState(np.array([1, 2, 3]), np.array([4, 5, 6]),
                       np.array([0., 1., 7.]), np.array([0., 1., 7.]),
                       np.array([0., 1., 7.]), np.array([0., 1., 7.]))

    # Test addition
    assert np.all((a + b).x == np.array([2, 4, 6]))
    assert np.all((a + b).y == np.array([8, 10, 12]))
    assert np.all((a + b).z == np.array([0., 2., 14.]))
    assert np.all((a + b).dx == np.array([0., 2., 14.]))
    assert np.all((a + b).dy == np.array([0., 2., 14.]))
    assert np.all((a + b).dz == np.array([0., 2., 14.]))

    # Test subtraction
    assert np.all((a - b).x == np.array([0, 0, 0]))
    assert np.all((a - b).y == np.array([0, 0, 0]))
    assert np.all((a - b).z == np.array([0., 0., 0.]))
    assert np.all((a - b).dx == np.array([0., 0., 0.]))
    assert np.all((a - b).dy == np.array([0., 0., 0.]))
    assert np.all((a - b).dz == np.array([0., 0., 0.]))

    # Addition of constant
    c = a + 1.
    assert np.all(c.x == np.array([2, 3, 4]))
    assert np.all(c.y == np.array([5, 6, 7]))
    assert np.all(c.z == np.array([1., 2., 8.]))
    assert np.all(c.dx == np.array([1., 2., 8.]))
    assert np.all(c.dy == np.array([1., 2., 8.]))
    assert np.all(c.dz == np.array([1., 2., 8.]))

    # Multiplication by constant
    c = a * 2.
    assert c is not None
    assert np.all(c.x == np.array([2, 4, 6]))
    assert np.all(c.y == np.array([8, 10, 12]))
    assert np.all(c.z == np.array([0., 2., 14.]))
    assert np.all(c.dx == np.array([0., 2., 14.]))
    assert np.all(c.dy == np.array([0., 2., 14.]))
    assert np.all(c.dz == np.array([0., 2., 14.]))

    # Multiplication by sequence
    c = a * [1, 2, 3, 4, 5, 6]
    assert c is not None
    assert np.all(c.x == np.array([1, 2, 3]))
    assert np.all(c.y == np.array([8, 10, 12]))
    assert np.all(c.z == np.array([0., 3., 21.]))
    assert np.all(c.dx == np.array([0., 4., 28.]))
    assert np.all(c.dy == np.array([0., 5., 35.]))
    assert np.all(c.dz == np.array([0., 6., 42.]))

    # Test append with arrays
    a.append(b)
    assert np.all(a.x == np.array([1, 2, 3, 1, 2, 3]))
    assert np.all(a.y == np.array([4, 5, 6, 4, 5, 6]))
    assert np.all(a.z == np.array([0., 1., 7., 0., 1., 7.]))
    assert np.all(a.dx == np.array([0., 1., 7., 0., 1., 7.]))
    assert np.all(a.dy == np.array([0., 1., 7., 0., 1., 7.]))
    assert np.all(a.dz == np.array([0., 1., 7., 0., 1., 7.]))

    # Test append with single elements
    c = CartesianState(1, 2, 3, 4, 5, 6)
    b.append(c)
    assert np.all(b.x == np.array([1, 2, 3, 1]))
    assert np.all(b.y == np.array([4, 5, 6, 2]))
    assert np.all(b.z == np.array([0., 1., 7., 3.]))
    assert np.all(b.dx == np.array([0., 1., 7., 4.]))
    assert np.all(b.dy == np.array([0., 1., 7., 5.]))
    assert np.all(b.dz == np.array([0., 1., 7., 6.]))

    return None


def test_KeplerianState_operations() -> None:

    a = KeplerianState(1, 2, 3, 0., 3.6, 0.9)
    b = KeplerianState(1, 2, 3, 0., 3.6, 0.9)

    # Test addition
    assert a + b == KeplerianState(2, 4, 6, 0., 7.2, 1.8)

    # Test subtraction
    assert a - b == KeplerianState(0, 0, 0, 0., 0., 0.)

    return None


def test_KeplerianStatesArray_operations() -> None:

    a = KeplerianState(np.array([1, 2, 3]), np.array([4, 5, 6]),
                       np.array([0., 1., 7.]), np.array([0., 1., 7.]),
                       np.array([0., 1., 7.]), np.array([0., 1., 7.]))

    b = KeplerianState(np.array([1, 2, 3]), np.array([4, 5, 6]),
                       np.array([0., 1., 7.]), np.array([0., 1., 7.]),
                       np.array([0., 1., 7.]), np.array([0., 1., 7.]))

    # Test addition
    assert np.all((a + b).a == np.array([2, 4, 6]))
    assert np.all((a + b).e == np.array([8, 10, 12]))
    assert np.all((a + b).i == np.array([0., 2., 14.]))
    assert np.all((a + b).Omega == np.array([0., 2., 14.]))
    assert np.all((a + b).omega == np.array([0., 2., 14.]))
    assert np.all((a + b).nu == np.array([0., 2., 14.]))

    return None


def test_SphericalGeocentricInit() -> None:

    # No corrections needed
    s = SphericalGeocentric(10., 0., 0.)
    assert s.r == 10.
    assert s.lat == 0.
    assert s.lon == 0.

    # Latitude over 90 deg
    s = SphericalGeocentric(10., 91., 0.)
    assert s.r == 10.
    assert s.lat == 89.
    assert s.lon == 180.

    s = SphericalGeocentric(10., 180., 0.)
    assert s.r == 10.
    assert s.lat == 0.
    assert s.lon == 180.

    s = SphericalGeocentric(10., 250., 0.)
    assert s.r == 10.
    assert s.lat == -70.
    assert s.lon == 180.

    # Latitude under -90 deg
    s = SphericalGeocentric(10., -91., 0.)
    assert s.r == 10.
    assert s.lat == -89.
    assert s.lon == 180.

    s = SphericalGeocentric(10., -180., 0.)
    assert s.r == 10.
    assert s.lat == 0.
    assert s.lon == 180.

    s = SphericalGeocentric(10., -250., 0.)
    assert s.r == 10.
    assert s.lat == 70.
    assert s.lon == 180.

    # Longitude over 180 deg
    s = SphericalGeocentric(10., 20., 181.)
    assert s.r == 10.
    assert s.lat == 20.
    assert s.lon == -179.

    s = SphericalGeocentric(10., 20., 270.)
    assert s.r == 10.
    assert s.lat == 20.
    assert s.lon == -90.

    s = SphericalGeocentric(10., 20., 360.)
    assert s.r == 10.
    assert s.lat == 20.
    assert s.lon == 0.

    # Longitude under -180 deg
    s = SphericalGeocentric(10., 20., -181.)
    assert s.r == 10.
    assert s.lat == 20.
    assert s.lon == 179.

    s = SphericalGeocentric(10., 20., -270.)
    assert s.r == 10.
    assert s.lat == 20.
    assert s.lon == 90.

    s = SphericalGeocentric(10., 20., -360.)
    assert s.r == 10.
    assert s.lat == 20.
    assert s.lon == 0.

    return None


def test_SphericalGeocentricSetters() -> None:

    s = SphericalGeocentric(10., 40., 60.)

    # Update values without errors
    s.set_r(20)
    assert s.r == 20.
    assert s.lat == 40.
    assert s.lon == 60.

    s.set_lat(190.)
    assert s.r == 20.
    assert s.lat == -10.
    assert s.lon == -120.

    s.set_lat(40.)
    assert s.r == 20.
    assert s.lat == 40.
    assert s.lon == -120.

    s.set_lon(60.)
    assert s.r == 20.
    assert s.lat == 40.
    assert s.lon == 60.


def test_SphericalGeocentricDifference() -> None:

    s1 = SphericalGeocentric(1., 90., 0.)
    s2 = SphericalGeocentric(1., -90., 0.)

    _, dlat, _ = s1.residual(s2)
    assert dlat == 180.

    s1.set_lon(50.)
    s2.set_lon(40.)
    assert s1.residual(s2)[2] == 10.
    assert s2.residual(s1)[2] == -10.

    s1.set_lon(190.)
    s2.set_lon(170.)
    assert s1.residual(s2)[2] == 20.
    assert s2.residual(s1)[2] == -20.

    s1.set_lon(370.)
    s2.set_lon(350.)
    assert s1.residual(s2)[2] == 20.
    assert s2.residual(s1)[2] == -20.

    return None
