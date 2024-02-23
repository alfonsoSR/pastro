from .base import BasePropagator
from ..types import CartesianState


class RungeKuttaPropagator(BasePropagator):

    def __init__(self, fun, t0, s0, tend, dt) -> None:

        super().__init__(fun, t0, s0, tend, dt)

        return None

    def _step_impl(self):

        # TODO: Handle units properly

        new_time_int, new_time_fr = self._update_epoch()
        h_int, h_fr = self._get_steps()

        k1 = self.fun(self.now_int, self.s, fr=self.now_fr)
        k2 = self.fun(self.now_int + 0.5 * self.h_int,
                      self.s + k1 * 0.5 * h_int + k1 * 0.5 * h_fr,
                      fr=self.now_fr + 0.5 * self.h_fr)
        k3 = self.fun(self.now_int + 0.5 * self.h_int,
                      self.s + k2 * 0.5 * h_int + k2 * 0.5 * h_fr,
                      fr=self.now_fr + 0.5 * self.h_fr)
        k4 = self.fun(self.now_int + self.h_int,
                      self.s + k3 * h_int + k3 * h_fr,
                      fr=self.now_fr + self.h_fr)
        sum_k = (k1 + k2 * 2. + k3 * 2. + k4) * (1. / 6.)

        new_s = self.s + sum_k * h_int + sum_k * h_fr
        assert isinstance(new_s, CartesianState)

        self.now_int = new_time_int
        self.now_fr = new_time_fr
        self.s = new_s

        return True
