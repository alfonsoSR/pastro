from ..types import CartesianState, Date
from ..utils import cglobal_diff
from typing import Type, overload
from ..solvers import BasePropagator
import numpy as np
from dataclasses import dataclass
from ..forces import ForceModel


@dataclass
class Individual:
    state: CartesianState
    loss: float


class Genetic:

    @overload
    def __init__(self, propagator: Type[BasePropagator],
                 forces: ForceModel, t0: Date, s0: CartesianState,
                 tend: Date, dt: float, s_ref: CartesianState,
                 seed: int) -> None:
        ...

    @overload
    def __init__(self, propagator: Type[BasePropagator],
                 forces: ForceModel, t0: float, s0: CartesianState,
                 tend: float, dt: float, s_ref: CartesianState,
                 seed: int) -> None:
        ...

    def __init__(self, propagator, forces, t0, s0, tend, dt, s_ref,
                 seed=859737) -> None:

        self.prop_class = propagator
        self.forces = forces
        self.t0 = t0
        self.s0 = s0
        self.tend = tend
        self.dt = dt
        self.ref = s_ref

        # Genetic algorithm parameters
        self.pool = np.linspace(-0.1, 0.1, 5)
        self.size = 10
        self.beta = 2

        # Random number generator
        self.gen = np.random.default_rng(seed)

        self.losslist = []

        return None

    def loss(self, s0: CartesianState) -> float:

        # Propagate orbit
        prop = self.prop_class(self.forces, self.t0, s0, self.tend, self.dt)
        _, _, s = prop.propagate()

        # Return cgdiff w.r.t reference
        return cglobal_diff(s, self.ref)[0]

    def mutate(self, population, prob):

        size = len(population)

        genes = self.gen.choice(self.pool, size=(size * 6)) + 1.
        options = np.array([genes, np.ones_like(genes)]).T

        mutate = np.reshape([self.gen.choice(opt, p=[prob, 1 - prob])
                             for opt in options], (size, 6))

        for idx, ind in enumerate(population):

            new_state = ind.state * mutate[idx]
            new_loss = self.loss(new_state)
            population[idx] = Individual(new_state, new_loss)

        return population

    def optimize(self, guess):

        population = [Individual(guess, 0) for _ in range(self.size)]

        population = self.mutate(population, 1)

        # Initialize offspring container
        offspring = [Individual(guess * 0, 1e30)
                     for _ in range(int(self.size * self.beta))]

        count = 0
        while True:

            # Sort population according to their loss
            population.sort(key=lambda x: x.loss)

            # Generate offspring
            for idx in range(0, self.size, 2):

                for jdx in range(int(2 * self.beta)):

                    # Recombination
                    _use = self.gen.choice([0, 1], size=(6))
                    use = np.array([_use, 1 - _use], dtype=float)

                    new_state = (population[idx].state * use[0] +
                                 population[idx + 1].state * use[1])
                    kdx = int(idx * self.beta + jdx)
                    offspring[kdx] = Individual(new_state, 0)

            # Mutation
            offspring = self.mutate(offspring, 0.05)

            # Combine parents and offspring
            all = population + offspring

            # Sort offspring according to their loss
            all.sort(key=lambda x: x.loss)

            # Replace population with offspring
            population = all[:self.size]

            self.losslist.append(population[0].loss)

            if count > 0:
                loss = population[0].loss
                old = self.losslist[count - 1]
                first = self.losslist[0]
                print(f"Step: {count} Loss: {loss:.2e}"
                      f" Rel: {(loss/old):.2f} First {(loss/first):.2e}")

            count += 1


def CGDLoss(s_ref: CartesianState, s_pred: CartesianState) -> float:
    """Loss function: Cumulative global difference

    :param s_ref: Reference state
    :param s_pred: Estimated state
    """

    ds = s_pred - s_ref
    gdiff = ds.r_mag

    return np.sqrt(np.sum(gdiff * gdiff, axis=0) / gdiff.shape[0])


def CGDLossGrad(s_ref: CartesianState, s_pred: CartesianState):
    """Gradient of loss function: Cumulative global difference"""

    param = 1. / (s_pred.x.shape[0] * CGDLoss(s_ref, s_pred))

    return (s_pred - s_ref) * param


class oldGradient:
    """Optimize initial state using gradient descent

    :param prop_cls: Propagator class
    :param setup: Propagator setup
    :param s_ref: Reference state
    """

    @overload
    def __init__(self, propagator: Type[BasePropagator],
                 forces: ForceModel, t0: Date, s0: CartesianState,
                 tend: Date, dt: float, s_ref: CartesianState,
                 alpha: float, pool_lim: float, mut_prob: float,
                 maxiter: int, seed: int) -> None:
        ...

    @overload
    def __init__(self, propagator: Type[BasePropagator],
                 forces: ForceModel, t0: float, s0: CartesianState,
                 tend: float, dt: float, s_ref: CartesianState,
                 alpha: float, pool_lim: float, mut_prob: float,
                 maxiter: int, seed: int) -> None:
        ...

    def __init__(self, propagator, forces, t0, s0, tend, dt, s_ref,
                 alpha=0.5, pool_lim=0.001, mut_prob=0.0001, maxiter=100,
                 seed=463672) -> None:

        # Orbit propagation
        self.prop_class = propagator
        self.forces = forces
        self.t0 = t0
        self.s0 = s0
        self.tend = tend
        self.dt = dt
        self.ref = s_ref

        # Parameters
        self.alpha = 50
        self.target = 10
        self.eps = 1e-3
        self.count = 0
        self.limit = 2
        self.max_iter = 2000

        # Caches
        self.last_loss = 1e20
        self.end = False
        self.best = CartesianState(0., 0., 0., 0., 0., 0.)
        self.loss_list = [0.0] * (self.max_iter + 1)

        # Random number generator
        self.gen = np.random.default_rng()

        return None

    def CGDLoss(self, s_pred: CartesianState) -> float:
        """Loss function: Cumulative global difference

        :param s_ref: Reference state
        :param s_pred: Estimated state
        """

        ds = s_pred - self.ref
        gdiff = ds.r_mag

        return np.sqrt(np.sum(gdiff * gdiff, axis=0) / gdiff.shape[0])

    def CGDLossGrad(self, s_pred: CartesianState) -> CartesianState:
        """Gradient of loss function: Cumulative global difference"""

        param: float = 1. / (s_pred.x.shape[0] * self.CGDLoss(s_pred))
        gradLoss = (s_pred - self.ref) * param
        assert isinstance(gradLoss, CartesianState)

        return gradLoss.initial_state

    def set_learning_rate(self, loss: float,
                          s0: CartesianState) -> CartesianState:
        """Adjust learning rate to accelerate convergence"""

        if loss > self.last_loss:
            if not (self.last_loss == 1e20):
                self.count += 1
                self.alpha *= 0.5
        else:
            self.best = s0

            if self.last_loss - loss < self.eps:
                delta = self.gen.uniform(-0.1, 0., (6,))
                s0 = CartesianState(s0.x + delta[0], s0.y + delta[1],
                                    s0.z + delta[2], s0.dx + delta[3],
                                    s0.dy + delta[4], s0.dz + delta[5])
            else:
                self.alpha = loss

        if self.count > self.limit:
            self.end = True

        return s0

    def optimize(self, s0: CartesianState) -> CartesianState:

        print("Optimizing initial state with gradient descent...")

        # Propagate orbit from initial guess
        prop = self.prop_class(self.forces, self.t0, s0, self.tend, self.dt)
        s = prop.propagate()[2]

        # Compute loss for initial guess
        loss = self.CGDLoss(s)

        count = 0
        self.alpha = loss
        self.loss_list[count] = loss

        while (not self.end) and (count < self.max_iter):

            # Adjust learning rate
            s0 = self.set_learning_rate(loss, s0)
            self.last_loss = loss

            grad = self.CGDLossGrad(s)
            s0 = s0 - grad * self.alpha
            prop = self.prop_class(self.forces, self.t0,
                                   s0, self.tend, self.dt)
            _, _, s = prop.propagate()
            loss = self.CGDLoss(s)
            count += 1
            self.loss_list[count] = loss

            print(f"Step: {count} Loss: {loss:.2f}")

        self.loss_list = self.loss_list[:count]

        return self.best
