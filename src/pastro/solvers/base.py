import numpy as np
from ..types import Date, CartesianState
from ..utils import date2jd
from ..forces import ForceModel
from typing import overload


class BasePropagator:

    @overload
    def __init__(self, fun: ForceModel, t0: Date, s0: CartesianState,
                 tend: Date, dt: float) -> None:
        ...

    @overload
    def __init__(self, fun: ForceModel, t0: float, s0: CartesianState,
                 tend: float, dt: float) -> None:
        ...

    def __init__(self, fun, t0, s0, tend, dt) -> None:

        # Evolution law
        self.fun = fun

        # Initial and final epochs and time step
        if isinstance(t0, Date) and isinstance(tend, Date):
            self.is_day = True
            self.now_int, self.now_fr = date2jd(t0, frac=True)
            self.tend_int, self.tend_fr = date2jd(tend, frac=True)
            dt /= 86400.
        elif isinstance(t0, float) and isinstance(tend, float):
            self.is_day = False
            self.now_int, self.now_fr = (t0, 0.)
            self.tend_int, self.tend_fr = (tend, 0.)
        else:
            raise TypeError("Unexpected type for initial or final epoch")

        self.now = self.now_int + self.now_fr
        self.tend = self.tend_int + self.tend_fr

        # Initial state
        if isinstance(s0, CartesianState):
            if not s0.is_scalar:
                raise TypeError(
                    "Components of initial state must be of scalar type")
        else:
            raise TypeError("Initial state must be of CartesianState type")

        self.s = s0

        # Time step
        delta_int = self.tend_int - self.now_int
        delta_fr = self.tend_fr - self.now_fr
        steps = (delta_int // dt) + (delta_fr // dt)
        self.h = (delta_int / steps) + (delta_fr / steps)
        self.h_int = np.trunc(self.h)
        self.h_fr = self.h - self.h_int

        if self.is_day:
            self.h_sec = self.h_int * 86400. + self.h_fr * 86400.
        else:
            self.h_sec = self.h

        # Termination criteria
        self.EPS = 1e-14
        self.done_int = max(self.EPS, 0.1 * self.h_int)
        self.done_fr = max(self.EPS, 0.1 * self.h_fr)

        # Status
        self.status = "running"

        return None

    def _update_epoch(self):

        new_fr = self.now_fr + self.h_fr
        int_part = np.trunc(new_fr)
        new_fr -= int_part
        new_int = self.now_int + self.h_int + int_part

        return new_int, new_fr

    def _get_steps(self):

        if self.is_day:
            return self.h_int * 86400, self.h_fr * 86400
        else:
            return self.h_int, self.h_fr

    def propagate(self, frac: bool = True):

        # Initialize output containers
        out_sol = CartesianState(self.s._x, self.s._y, self.s._z,
                                 self.s._dx, self.s._dy, self.s._dz)
        out_time = [self.now_int]
        out_time_fr = [self.now_fr]

        # Main loop
        prop_status = None
        while prop_status is None:

            # Perform propagation step
            self.step()

            # Check status after step
            if self.status == "finished":
                prop_status = 0
            elif self.status == "failed":
                prop_status = -1
                break

            # Update output arrays
            out_sol.append(self.s)
            out_time.append(self.now_int)
            out_time_fr.append(self.now_fr)

        # Generate output
        time_int = np.array(out_time, dtype=np.float64)
        time_fr = np.array(out_time_fr, dtype=np.float64)
        time = np.array(time_int + time_fr, dtype=np.float64)

        if frac:
            return time_int, time_fr, out_sol
        else:
            return time, 0. * time_fr, out_sol

    def step(self):

        # Check status before performing step
        if self.status != "running":
            raise RuntimeError("Attempted step on non-running propagator")

        # Perform step
        success = self._step_impl()

        # Check for termination conditions
        if not success:
            self.status = "failed"
        else:
            if ((self.tend_int - self.now_int < self.done_int) and
                    (self.tend_fr - self.now_fr < self.done_fr)):
                self.status = "finished"

        return None

    def _step_impl(self):
        raise NotImplementedError
