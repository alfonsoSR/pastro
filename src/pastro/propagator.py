import requests
from sgp4.api import Satrec
from .utils import date2jd
from .types import Date, Vector, CartesianState
import numpy as np
from typing import overload


class SGP4Propagator:
    """Propagate satellite orbit from its TLE

    :param norad_id: NORAD ID of the satellite
    :param tle: TLE as list of strings (Name, line 1 and line 2) or None to
    retrieve latest from Celestrak (Default: None).
    """

    def __init__(self, norad_id: int, tle: list[str] | None = None) -> None:

        self.id = norad_id

        # Read provided TLE or retrieve from Celestrak
        if tle is None:
            self.tle = self._get_tle()
        else:
            self.tle = tle

        # Initialize satellite object with TLE
        self.sat = Satrec.twoline2rv(self.tle[1], self.tle[2])

        return None

    def _get_tle(self) -> list[str]:
        """Retrieve latest TLE from Celestrak for selected satellite

        :return: TLE as list of three strings (Name, line 1 and line 2)
        """

        # Retrieve TLE from Celestrak
        tle = requests.get(
            f"https://celestrak.org/NORAD/elements/gp.php?CATNR={self.id}",
            timeout=5, verify=True)
        tle.raise_for_status()

        # Ensure that we got TLE data
        if tle.text == "No GP data found":
            raise ValueError("No TLE data found for given NORAD ID")

        return tle.text.splitlines()

    def relative_time(self, date: Date,
                      frac: bool = False) -> tuple[float, float]:
        """Difference between given epoch and TLE of propagator in seconds.

        :param date: Epoch to compare with TLE.
        :param frac: Whether to return the integer and fractional part of the
        difference separately. If set to False, the whole difference is
        returned in the first element of the tuple and the second is set to
        zero (default: False).
        """

        # Get TLE epoch as Julian date
        tle_jd = self.sat.jdsatepoch
        tle_jdf = self.sat.jdsatepochF

        # Get given epoch as Julian date
        date_jd, date_jdf = date2jd(date, frac=True)

        # Calculate difference in seconds
        d_jd = (date_jd - tle_jd) * 86400.
        d_jdf = (date_jdf - tle_jdf) * 86400.

        # Return based on selected format
        if frac:
            return d_jd, d_jdf
        else:
            return d_jd + d_jdf, 0.

    @overload
    def _generate_epochs(self, t0: float, tend: float,
                         step: float) -> tuple[Vector, Vector]:
        ...

    @overload
    def _generate_epochs(self, t0: Date, tend: Date,
                         step: float) -> tuple[Vector, Vector]:
        ...

    def _generate_epochs(self, t0, tend, step):
        """Generate range of epochs from initial and final epochs and step

        If the epochs are passed as Date objects, the returned tuple contains
        the integer and fractional parts of the Julian date epochs separately.
        Otherwise, the first element of the tuple contains the epochs in
        seconds past t0 and the second element is empty.

        :param t0: Initial epoch.
        :param tend: Final epoch.
        :param step: Time step in seconds.
        """

        if t0 == tend:
            return date2jd(t0, frac=True)

        if isinstance(t0, Date) and isinstance(tend, Date):

            jd0, fr0 = date2jd(t0, frac=True)
            jd_end, fr_end = date2jd(tend, frac=True)
            step /= 86400.

            # Find number of steps with original step
            delta_jd = (jd_end - jd0)
            delta_fr = (fr_end - fr0)
            steps = (delta_jd // step) + (delta_fr // step)
            isteps = int(steps)

            step_jd = delta_jd / steps
            step_fr = delta_fr / steps
            step = step_jd + step_fr

            time_jd = np.array([jd0 + i * step_jd for i in range(isteps + 1)])
            time_fr = np.array([fr0 + i * step_fr for i in range(isteps + 1)])

            return time_jd, time_fr

        elif isinstance(t0, float) and isinstance(tend, float):

            steps = int((tend - t0) // step)
            time = np.linspace(t0, tend, steps)

            return time, 0. * time

    def propagate(self, t0: Date, tend: Date,
                  dt: float) -> tuple[Vector, Vector, CartesianState]:
        """Propagate orbit from TLE between given epochs using SGP4.

        :param t0: Initial epoch.
        :param tend: Final epoch.
        :param dt: Time step in seconds.
        :return time_jd: Integer part of epochs in which the state is known.
        :return time_fr: Fractional part of epochs in which the state is known.
        :return s: Cartesian state vector [m, m/s]
        """

        time_jd, time_fr = self._generate_epochs(t0, tend, dt)
        e, r, u = self.sat.sgp4_array(time_jd, time_fr)

        # Ensure that propagation succeeded
        if e.any():
            raise ValueError("Propagation failed")

        # Return Cartesian state vectors
        r *= 1e3
        u *= 1e3
        s = CartesianState(r[:, 0], r[:, 1], r[:, 2],
                           u[:, 0], u[:, 1], u[:, 2])

        return time_jd, time_fr, s
