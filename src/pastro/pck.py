from dataclasses import dataclass


@dataclass
class Earth:

    R: float = 6.371e6
    mu: float = 3.986004418e14
    omega: float = 71.92115e-6
    j2: float = 1.083e-3


@dataclass
class Moon:

    mu: float = 4.9018000066e12


@dataclass
class Sun:

    R: float = 696e6
    mu: float = 1.32712440018e20
