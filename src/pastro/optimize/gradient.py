from typing import Type, overload
from ..solvers import BasePropagator
from ..forces import ForceModel
import numpy as np
from ..types import CartesianState, Date

# class NewGradient:
    
#     @overload
#     def __init__(self, propagator: Type[BasePropagator],
#                  forces: ForceModel, t0: Date, s0: CartesianState,
#                  tend: Date, dt: float, s_ref: CartesianState,
#                  alpha: float, pool_lim: float, mut_prob: float,
#                  maxiter: int, seed: int) -> None:
#         ...

#     @overload
#     def __init__(self, propagator: Type[BasePropagator],
#                  forces: ForceModel, t0: float, s0: CartesianState,
#                  tend: float, dt: float, s_ref: CartesianState,
#                  alpha: float, pool_lim: float, mut_prob: float,
#                  maxiter: int, seed: int) -> None:
#         ...

#     def __init__(self, propagator, forces, t0, s0, tend, dt, s_ref,
#                  alpha=0.85, pool_lim=0.01, mut_prob=0.0001, maxiter=100,
#                  seed=463672) -> None:
        
#         # Propagator setup
#         self.type = propagator
#         self.forces = forces
#         self.t0 = t0
#         self.tend = tend
#         self.s0 = s0
#         self.dt = dt
#         self.ref = s_ref			# Reference orbit: SGP4
        
#         # Newton method parameters
#         self.alpha = alpha
        
#         return None
    
#     def cgd_loss(self, s: CartesianState):
        
#         ds = s - self.ref
#         gdiff = ds.r_mag
        
#         return np.sqrt(np.sum(gdiff * gdiff) / gdiff.shape[0])
    
#     def cgd_loss_grad(self, s: CartesianState):
        
#         ds = (s - self.ref).asarray()
        
#         ds = (s - self.ref).asarray().swapaxes(0, 1)
#         return ds / (s.x.shape[0] * self.cgd_loss(s))
    
#     def cgd_loss_hessian(self, s: CartesianState):
        
#         loss = self.cgd_loss(s)
#         nL = s.x.shape[0] * loss
#         nL2 = nL * loss
#         inv_nL = 1. / nL
#         inv_nL2 = 1. / nL2
#         ds = (s - self.ref).asarray().swapaxes(0, 1)
        
#         out = inv_nL * (np.identity(6)[None, :, :] -
#                          ds[:, :, None] * ds[:, None, :] * inv_nL2)
#         return out
    
#     def optimize(self, guess: CartesianState):
        
#         # Calculate loss for initial guess
#         prop = self.type(self.forces, self.t0, guess, self.tend, self.dt)
#         s = prop.propagate()[2]
#         loss = self.cgd_loss(s)
#         print(f"Loss: {loss:.2f}")
        
#         # Optimization loop
#         run = True
#         count = 0
#         while run:
            
#             if count > 100:
#                 run = False
#                 continue
            
#             # Update guess
#             self.cgd_loss_grad(s)
#             self.cgd_loss_hessian(s)
#             exit(0)
            
            


class Gradient:

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
                 alpha=0.85, pool_lim=0.01, mut_prob=0.0001, maxiter=100,
                 seed=463672) -> None:

        # Propagator setup
        self.type = propagator
        self.forces = forces
        self.t0 = t0
        self.tend = tend
        self.s0 = s0
        self.dt = dt
        self.ref = s_ref

        # Optimization parameters
        self.alpha = alpha
        self.mut_prob = mut_prob

        # Termination conditions
        self.run = True
        self.eps = 0.1
        self.close_limit = 1
        self.maxiter = maxiter

        # Counters and switches
        self.force = 0
        self.count_close = 0

        # Caches
        self.ref_s0 = self.ref.initial_state
        self.old_loss = 0.

        ref_prop = self.type(self.forces, self.t0, self.s0, self.tend, self.dt)
        s_ref_prop = ref_prop.propagate()[2]
        self.ref_loss = self.CGDLoss(s_ref_prop)
        
        self.loss_cache = np.zeros((self.maxiter + 1, 1), dtype=float)

        # Misc
        self.rand = np.random.default_rng(seed)
        self.pool = np.linspace(-pool_lim, pool_lim, 5)

        return None
    
    def loss_from_perturbation(self, ds: CartesianState) -> float:
        """Calculate loss from perturbation of initial state"""
        
        return np.sqrt(
            np.sum(ds.r_mag * ds.r_mag, axis=0) / ds.r_mag.shape[0])

    def CGDLoss(self, s: CartesianState) -> float:
        """Calculate loss with respect to reference state"""

        ds = s - self.ref
        gdiff = ds.r_mag

        return np.sqrt(np.sum(gdiff * gdiff, axis=0) / gdiff.shape[0])

    def CGDLossGrad(self, s: CartesianState) -> np.ndarray:
        """Calculate gradient of loss function at given state"""

        # ds = s - self.ref
        ds = (s.initial_state - self.ref_s0).asarray()
        param = 1. / (s.x.shape[0] * self.CGDLoss(s))
        grad = ds * param

        return grad

        # return grad.asarray().T

    def CGDLossHessian(self, s: CartesianState) -> np.ndarray:

        loss = self.CGDLoss(s)
        nL = s.x.shape[0] * loss
        nL2 = nL * loss
        inv_nL = 1. / nL
        inv_nL2 = 1. / nL2
        ds = (s.initial_state - self.ref_s0).asarray()

        return inv_nL * (np.identity(6) - ds[:, None] * ds[None, :] * inv_nL2)

    def mutate(self, s: CartesianState):

        params = self.rand.choice(self.pool, 6) + 1.

        if self.force:
            out = s * params
            self.force = 0
        else:
            base = np.ones_like(params)

            use = self.rand.choice([0, 1], 6,
                                   p=[self.mut_prob, 1 - self.mut_prob])

            out = s * (base * use + params * (1 - use))

        assert isinstance(out, CartesianState)
        return out

    def step_size(self, s0: CartesianState, loss: float):
        """Select step size based on current guess and loss"""

        return None

    def optimize(self, guess: CartesianState):

        # Calculate loss for initial guess
        prop = self.type(self.forces, self.t0, guess, self.tend, self.dt)
        s = prop.propagate()[2]

        # Compute loss for initial guess
        loss = self.CGDLoss(s)
        print(f"Loss: {loss:.2f}")

        # Optimization loop
        count = 0
        while self.run:

            if count > self.maxiter:
                self.run = False
                continue

            self.old_loss = loss

            if self.count_close == self.close_limit:
                self.count_close = 0
                self.force = 1

            # Update guess
            h = - np.linalg.solve(self.CGDLossHessian(s), self.CGDLossGrad(s))
            guess = guess + h * self.alpha
            guess = self.mutate(guess)

            # Calculate loss for updated guess
            prop = self.type(self.forces, self.t0, guess, self.tend, self.dt)
            s = prop.propagate()[2]
            loss = self.CGDLoss(s)
            self.loss_cache[count] = loss

            # Update target loss if we improve original target
            if loss < self.ref_loss:
                self.ref_loss = loss
                self.ref_s0 = guess

            # Mutate initial state if method converges
            if np.abs(self.old_loss - loss) < self.eps:
                self.count_close += 1

            count += 1

            print(f"Step: {count} Loss: {loss:.2f} Count: {self.count_close}")

        # Select step-size

        # Perform optimization
        
        print("End of optimization")
        print(f"Best loss: {self.ref_loss:.2f}")

        return self.ref_s0
