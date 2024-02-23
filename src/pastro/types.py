from dataclasses import dataclass
from numpy.typing import NDArray
import numpy as np
from typing import Sequence, overload

# Type alias for arrays of floats
Vector = NDArray[np.float64]
Scalar = int | float


# Date type
@dataclass
class Date:
    """Date and time in Gregorian calendar.

    :param day:
    :param month:
    :param year:
    :param hour:
    :param minute:
    :param second:
    """
    day: int
    month: int
    year: int
    hour: int
    minute: int
    second: int

    def tolist(self) -> list[int]:
        """Return date as list of integers.

        NOTE: Default format is [year, month, day, hour, minute, second]
        """
        return [self.year, self.month, self.day,
                self.hour, self.minute, self.second]


class CartesianState:

    def __init__(self, x, y, z, dx, dy, dz) -> None:

        if isinstance(x, Scalar):
            self.is_scalar = True

            assert isinstance(y, Scalar)
            assert isinstance(z, Scalar)
            assert isinstance(dx, Scalar)
            assert isinstance(dy, Scalar)
            assert isinstance(dz, Scalar)

            self._x = np.array([x], dtype=np.float64)
            self._y = np.array([y], dtype=np.float64)
            self._z = np.array([z], dtype=np.float64)
            self._dx = np.array([dx], dtype=np.float64)
            self._dy = np.array([dy], dtype=np.float64)
            self._dz = np.array([dz], dtype=np.float64)

        elif isinstance(x, np.ndarray):
            self.is_scalar = False

            assert isinstance(y, np.ndarray)
            assert isinstance(z, np.ndarray)
            assert isinstance(dx, np.ndarray)
            assert isinstance(dy, np.ndarray)
            assert isinstance(dz, np.ndarray)

            assert x.size == y.size
            assert x.size == z.size
            assert x.size == dx.size
            assert x.size == dy.size
            assert x.size == dz.size

            self._x = x
            self._y = y
            self._z = z
            self._dx = dx
            self._dy = dy
            self._dz = dz

        else:
            raise TypeError(
                "Arguments must be scalars or vectors of equal size")

        return None

    @property
    def x(self):
        return self._x[0] if self.is_scalar else self._x

    @property
    def y(self):
        return self._y[0] if self.is_scalar else self._y

    @property
    def z(self):
        return self._z[0] if self.is_scalar else self._z

    @property
    def dx(self):
        return self._dx[0] if self.is_scalar else self._dx

    @property
    def dy(self):
        return self._dy[0] if self.is_scalar else self._dy

    @property
    def dz(self):
        return self._dz[0] if self.is_scalar else self._dz

    @property
    def r(self):
        """Position vector as numpy array."""
        return np.array([self.x, self.y, self.z])

    @property
    def u(self):
        """Velocity vector as numpy array."""
        return np.array([self.dx, self.dy, self.dz])

    @property
    def r_mag(self):
        """Magnitude of position vector."""
        return np.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    @property
    def u_mag(self):
        """Magnitude of velocity vector."""
        return np.sqrt(self.dx * self.dx
                       + self.dy * self.dy + self.dz * self.dz)

    def __add__(self, other):

        if isinstance(other, CartesianState):
            return CartesianState(self.x + other.x, self.y + other.y,
                                  self.z + other.z, self.dx + other.dx,
                                  self.dy + other.dy, self.dz + other.dz)

        elif isinstance(other, int | float):
            return CartesianState(self.x + other, self.y + other,
                                  self.z + other, self.dx + other,
                                  self.dy + other, self.dz + other)

        elif isinstance(other, CartesianStateDerivative):
            # Assume that the derivative has been multiplied by some constant
            # so that the physics make sense
            return CartesianState(self.x + other.dx, self.y + other.dy,
                                  self.z + other.dz, self.dx + other.ddx,
                                  self.dy + other.ddy, self.dz + other.ddz)

        elif (isinstance(other, np.ndarray) and (self.is_scalar)
              and (other.shape == (6,))):
            return CartesianState(self.x + other[0], self.y + other[1],
                                  self.z + other[2], self.dx + other[3],
                                  self.dy + other[4], self.dz + other[5])
        else:
            raise TypeError(
                f"Cannot add object of type {type(other).__name__}")

    def __sub__(self, other):

        if isinstance(other, CartesianState):
            return CartesianState(self.x - other.x, self.y - other.y,
                                  self.z - other.z, self.dx - other.dx,
                                  self.dy - other.dy, self.dz - other.dz)

        elif isinstance(other, int | float):
            return CartesianState(self.x - other, self.y - other,
                                  self.z - other, self.dx - other,
                                  self.dy - other, self.dz - other)

        elif isinstance(other, CartesianStateDerivative):
            return CartesianState(self.x - other.dx, self.y - other.dy,
                                  self.z - other.dz, self.dx - other.ddx,
                                  self.dy - other.ddy, self.dz - other.ddz)

        else:
            raise TypeError(
                f"Cannot subtract object of type {type(other).__name__}")

    def __mul__(self, other):

        if isinstance(other, int | float):
            return CartesianState(self.x * other, self.y * other,
                                  self.z * other, self.dx * other,
                                  self.dy * other, self.dz * other)

        elif isinstance(other, Sequence | np.ndarray):
            if isinstance(other[0], int | float) and (len(other) == 6):
                return CartesianState(self.x * other[0], self.y * other[1],
                                      self.z * other[2], self.dx * other[3],
                                      self.dy * other[4], self.dz * other[5])
        else:
            raise TypeError(
                f"Cannot multiply object of type {type(other).__name__}")

    def __eq__(self, other):
        return ((self.x == other.x) and (self.y == other.y) and
                (self.z == other.z) and (self.dx == other.dx) and
                (self.dy == other.dy) and (self.dz == other.dz))

    def __iter__(self):

        return iter([CartesianState(self._x[i], self._y[i], self._z[i],
                                    self._dx[i], self._dy[i], self._dz[i])
                     for i in range(len(self._x))])

        # return iter([self.x, self.y, self.z, self.dx, self.dy, self.dz])

    def __repr__(self) -> str:
        return f"""
    CartesianState(x = {self.x},
                   y = {self.y},
                   z = {self.z},
                   dx = {self.dx},
                   dy = {self.dy},
                   dz = {self.dz})"""

    def append(self, other) -> None:
        """Append other state to current state."""
        if self.is_scalar:
            self.is_scalar = False
            # raise TypeError("Cannot append to scalar state")

        self._x = np.append(self._x, other._x)
        self._y = np.append(self._y, other._y)
        self._z = np.append(self._z, other._z)
        self._dx = np.append(self._dx, other._dx)
        self._dy = np.append(self._dy, other._dy)
        self._dz = np.append(self._dz, other._dz)

        return None

    def asarray(self):
        """Return state as numpy array"""
        return np.array([self.x, self.y, self.z, self.dx, self.dy, self.dz])

    @property
    def initial_state(self):
        """Return initial state"""
        return CartesianState(self._x[0], self._y[0], self._z[0],
                              self._dx[0], self._dy[0], self._dz[0])


class CartesianStateDerivative:

    def __init__(self, dx, dy, dz, ddx, ddy, ddz) -> None:

        if isinstance(dx, Scalar):
            self.is_scalar = True

            assert isinstance(dy, Scalar)
            assert isinstance(dz, Scalar)
            assert isinstance(ddx, Scalar)
            assert isinstance(ddy, Scalar)
            assert isinstance(ddz, Scalar)

            self._dx = np.array([dx], dtype=np.float64)
            self._dy = np.array([dy], dtype=np.float64)
            self._dz = np.array([dz], dtype=np.float64)
            self._ddx = np.array([ddx], dtype=np.float64)
            self._ddy = np.array([ddy], dtype=np.float64)
            self._ddz = np.array([ddz], dtype=np.float64)

        elif isinstance(dx, np.ndarray):
            self.is_scalar = False

            assert isinstance(dy, np.ndarray)
            assert isinstance(dz, np.ndarray)
            assert isinstance(ddx, np.ndarray)
            assert isinstance(ddy, np.ndarray)
            assert isinstance(ddz, np.ndarray)

            assert dx.size == dy.size
            assert dx.size == dz.size
            assert dx.size == ddx.size
            assert dx.size == ddy.size
            assert dx.size == ddz.size

            self._dx = dx
            self._dy = dy
            self._dz = dz
            self._ddx = ddx
            self._ddy = ddy
            self._ddz = ddz

        else:
            raise TypeError(
                "Arguments must be scalars or vectors of equal size")

        return None

    @property
    def dx(self):
        return self._dx[0] if self.is_scalar else self._dx

    @property
    def dy(self):
        return self._dy[0] if self.is_scalar else self._dy

    @property
    def dz(self):
        return self._dz[0] if self.is_scalar else self._dz

    @property
    def ddx(self):
        return self._ddx[0] if self.is_scalar else self._ddx

    @property
    def ddy(self):
        return self._ddy[0] if self.is_scalar else self._ddy

    @property
    def ddz(self):
        return self._ddz[0] if self.is_scalar else self._ddz

    @property
    def u(self):
        """Velocity vector as numpy array."""
        return np.array([self.dx, self.dy, self.dz])

    @property
    def a(self):
        """Acceleration vector as numpy array."""
        return np.array([self.ddx, self.ddy, self.ddz])

    @property
    def u_mag(self):
        """Magnitude of velocity vector."""
        return np.sqrt(self.dx * self.dx
                       + self.dy * self.dy + self.dz * self.dz)

    @property
    def a_mag(self):
        """Magnitude of acceleration vector."""
        return np.sqrt(self.ddx * self.ddx
                       + self.ddy * self.ddy + self.ddz * self.ddz)

    def __add__(self, other):

        if isinstance(other, CartesianStateDerivative):
            return CartesianStateDerivative(
                self.dx + other.dx, self.dy + other.dy,
                self.dz + other.dz, self.ddx + other.ddx,
                self.ddy + other.ddy, self.ddz + other.ddz)

        elif isinstance(other, int | float):
            return CartesianStateDerivative(
                self.dx + other, self.dy + other,
                self.dz + other, self.ddx + other,
                self.ddy + other, self.ddz + other)

        elif isinstance(other, CartesianState):
            # Assume the cartesian state is multiplied by some constant so
            # that the physics make sense
            return CartesianStateDerivative(
                self.dx + other.dx, self.dy + other.dy,
                self.dz + other.dz, self.ddx + other.dx,
                self.ddy + other.dy, self.ddz + other.dz)

        else:
            raise TypeError(
                f"Cannot add object of type {type(other).__name__}")

    def __sub__(self, other):
        return CartesianStateDerivative(
            self.dx - other.dx, self.dy - other.dy, self.dz - other.dz,
            self.ddx - other.ddx, self.ddy - other.ddy, self.ddz - other.ddz)

    def __mul__(self, other):

        if isinstance(other, int | float):
            return CartesianStateDerivative(self.dx * other, self.dy * other,
                                            self.dz * other, self.ddx * other,
                                            self.ddy * other, self.ddz * other)
        elif isinstance(other, Sequence):
            if isinstance(other[0], int | float) and (len(other) == 6):
                return CartesianStateDerivative(
                    self.dx * other[0], self.dy * other[1],
                    self.dz * other[2], self.ddx * other[3],
                    self.ddy * other[4], self.ddz * other[5])
        else:
            raise TypeError(
                f"Cannot multiply object of type {type(other).__name__}")

    def __eq__(self, other):
        return ((self.dx == other.dx) and (self.dy == other.dy) and
                (self.dz == other.dz) and (self.ddx == other.ddx) and
                (self.ddy == other.ddy) and (self.ddz == other.ddz))

    def __iter__(self):
        return iter([self.dx, self.dy, self.dz, self.ddx, self.ddy, self.ddz])

    def append(self, other) -> None:
        """Append other state to current state."""
        if self.is_scalar:
            self.is_scalar = False
            # raise TypeError("Cannot append to scalar state")

        self._dx = np.append(self._dx, other._dx)
        self._dy = np.append(self._dy, other._dy)
        self._dz = np.append(self._dz, other._dz)
        self._ddx = np.append(self._ddx, other._ddx)
        self._ddy = np.append(self._ddy, other._ddy)
        self._ddz = np.append(self._ddz, other._ddz)

        return None

    def add_acceleration(self, other) -> None:
        """Add accelerations of two state vector derivatives

        Updates the accelerations by adding those of the other state derivative
        and preserves the original velocities.

        Solves the problem of combining the contributions of multiple forces.
        """

        if not isinstance(other, CartesianStateDerivative):
            raise TypeError("Input must be of type CartesianStateDerivative")

        self._ddx += other._ddx
        self._ddy += other._ddy
        self._ddz += other._ddz

        return None

    @property
    def initial_state(self):
        """Return initial state"""
        return CartesianStateDerivative(
            self._dx[0], self._dy[0], self._dz[0],
            self._ddx[0], self._ddy[0], self._ddz[0])


class KeplerianState:
    """State vector in keplerian orbital elements.

    :param a: Semi-major axis or time series of them [m]
    :param e: Eccentricity or time series of them [-]
    :param i: Inclination or time series of them [deg]
    :param Omega: Longitude of the ascending node or time series of them [deg]
    :param omega: Argument of periapsis or time series of them [deg]
    :param nu: True anomaly or time series of them [deg]
    """

    def __init__(self, a: Scalar | Vector, e: Scalar | Vector,
                 i: Scalar | Vector, Omega: Scalar | Vector,
                 omega: Scalar | Vector, nu: Scalar | Vector) -> None:
        """Turn input into arrays and ensure they all have the same size."""

        if isinstance(a, Scalar):
            a = np.array([a])
        if isinstance(e, Scalar):
            e = np.array([e])
        if isinstance(i, Scalar):
            i = np.array([i])
        if isinstance(Omega, Scalar):
            Omega = np.array([Omega])
        if isinstance(omega, Scalar):
            omega = np.array([omega])
        if isinstance(nu, Scalar):
            nu = np.array([nu])

        ref_size = a.size
        assert e.size == ref_size
        assert i.size == ref_size
        assert Omega.size == ref_size
        assert omega.size == ref_size
        assert nu.size == ref_size

        # Ensure that all angles are between 0 and 360 degrees
        for idx, _ in enumerate(i):
            while i[idx] > 360.:
                i[idx] -= 360.
            while i[idx] < 0.:
                i[idx] += 360.

        for idx, _ in enumerate(Omega):
            while Omega[idx] > 360.:
                Omega[idx] -= 360
            while Omega[idx] < 0.:
                Omega[idx] += 360

        for idx, _ in enumerate(omega):
            while omega[idx] > 360.:
                omega[idx] -= 360.
            while omega[idx] < 0.:
                omega[idx] += 360.

        for idx, _ in enumerate(nu):
            while nu[idx] > 360.:
                nu[idx] -= 360,
            while nu[idx] < 0.:
                nu[idx] += 360.

        self.a = a
        self.e = e
        self.i = i
        self.Omega = Omega
        self.omega = omega
        self.nu = nu

        return None

    def __add__(self, other):
        print("WARNING: Addition of angles is not properly implemented")

        return KeplerianState(self.a + other.a, self.e + other.e,
                              self.i + other.i, self.Omega + other.Omega,
                              self.omega + other.omega, self.nu + other.nu)

    def __sub__(self, other):
        print("WARNING: Subtraction of angles is not properly implemented")

        return KeplerianState(self.a - other.a, self.e - other.e,
                              self.i - other.i, self.Omega - other.Omega,
                              self.omega - other.omega, self.nu - other.nu)

    def __eq__(self, other):
        return ((self.a == other.a) and (self.e == other.e) and
                (self.i == other.i) and (self.Omega == other.Omega) and
                (self.omega == other.omega) and (self.nu == other.nu))


class SphericalCoordinates:
    """Position vector in spherical coordinates

    NOTE: All angles are converted to radians during initialization.

    :param r: Radius or time series of them [m]
    :param lat: Latitude or time series of them [rad]
    :param long: Longitude or time series of them [rad]
    """

    def __init__(self, r: Scalar | Vector, lat: Scalar | Vector,
                 long: Scalar | Vector) -> None:
        """Turn input into arrays and ensure they all have the same size."""

        if isinstance(r, Scalar):
            r = np.array([r])
        if isinstance(lat, Scalar):
            lat = np.array([lat])
        if isinstance(long, Scalar):
            long = np.array([long])

        ref_size = r.size
        assert lat.size == ref_size
        assert long.size == ref_size

        # Ensure that all angles are in range
        for idx, _ in enumerate(lat):
            while lat[idx] > 90.:
                lat[idx] -= 90.
            while lat[idx] < -90.:
                lat[idx] += 90.

        for idx, _ in enumerate(long):
            while long[idx] > 180.:
                long[idx] -= 180.
            while long[idx] < -180:
                long[idx] += 180.

        self.r = r
        self.lat = lat
        self.long = long

        return None


class SphericalGeocentric:
    """Position vector in geocentric radius, latitude and longitude

    :param r: Orbital radius [m]
    :param lat: Geocentric latitude [deg]
    :param lon: Geocentric longitude [deg]
    """

    @overload
    def __init__(self, r: Scalar, lat: Scalar, lon: Scalar) -> None: ...

    @overload
    def __init__(self, r: Vector, lat: Vector, lon: Vector) -> None: ...

    def __init__(self, r, lat, lon) -> None:

        if isinstance(r, Scalar):
            self.is_scalar = True

            assert isinstance(lat, Scalar)
            assert isinstance(lon, Scalar)

            self._r = np.array([r], dtype=float)
            self._lat = np.array([lat], dtype=float)
            self._lon = np.array([lon], dtype=float)

        elif isinstance(r, np.ndarray):
            self.is_scalar = False

            assert isinstance(lat, np.ndarray)
            assert isinstance(lon, np.ndarray)

            assert r.size == lat.size
            assert r.size == lon.size

            self._r = r
            self._lat = lat
            self._lon = lon

        else:
            raise TypeError(
                "Arguments must be scalar or vectors of equal size")

        # Ensure that all angles are within range
        self._constraint_angles()

        return None

    @property
    def r(self):
        return self._r[0] if self.is_scalar else self._r

    @property
    def lat(self):
        return self._lat[0] if self.is_scalar else self._lat

    @property
    def lon(self):
        return self._lon[0] if self.is_scalar else self._lon

    @overload
    def set_r(self, r: Scalar) -> None: ...

    @overload
    def set_r(self, r: Vector) -> None: ...

    def set_r(self, r) -> None:

        if isinstance(r, Scalar) and self.is_scalar:
            self._r[0] = r
        elif isinstance(r, np.ndarray) and (not self.is_scalar):
            if self._r.size == r.size:
                self._r = r
            else:
                raise ValueError(
                    "New and old radii sequences must have the same size")
        else:
            raise TypeError("New and old radius must be of same type")

        return None

    @overload
    def set_lat(self, lat: Scalar) -> None: ...

    @overload
    def set_lat(self, lat: Vector) -> None: ...

    def set_lat(self, lat) -> None:

        if isinstance(lat, Scalar) and self.is_scalar:
            self._lat[0] = lat
        elif isinstance(lat, np.ndarray) and (not self.is_scalar):
            if self._lat.size == lat.size:
                self._lat = lat
            else:
                raise ValueError(
                    "New and old latitude sequences must have the same size")
        else:
            raise TypeError("New and old latitudes must be of same type")

        self._constraint_angles()

        return None

    @overload
    def set_lon(self, lon: Scalar) -> None: ...

    @overload
    def set_lon(self, lon: Vector) -> None: ...

    def set_lon(self, lon) -> None:

        if isinstance(lon, Scalar) and self.is_scalar:
            self._lon[0] = lon
        elif isinstance(lon, np.ndarray) and (not self.is_scalar):
            if self._lon.size == lon.size:
                self._lon = lon
            else:
                raise ValueError(
                    "New and old longitude sequences must have the same size")
        else:
            raise TypeError("New and old longitudes must be of same type")

        self._constraint_angles()

        return None

    def _constraint_angles(self) -> None:

        for idx, _ in enumerate(self._lat):

            # Correct latitude and flip longitude if needed
            while (self._lat[idx] > 90.) or (self._lat[idx] < -90.):
                if self._lat[idx] > 90.:
                    self._lat[idx] = 180. - self._lat[idx]
                elif self._lat[idx] < -90.:
                    self._lat[idx] = -180. - self._lat[idx]

                self._lon[idx] += 180.

            # Correct longitude
            while (self._lon[idx] > 180.) or (self._lon[idx] < -180.):
                if self._lon[idx] > 180.:
                    self._lon[idx] -= 360.
                elif self._lon[idx] < -180.:
                    self._lon[idx] += 360.

        return None

    def residual(self, other):

        dr = self.r - other.r
        dlat = self.lat - other.lat

        lons = np.unwrap([self.lon, other.lon], period=360.)
        dlon = lons[0] - lons[1]

        return dr, dlat, dlon
