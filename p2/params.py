# All inputs live here. Tags: [DESIGN] from CAD, not measured yet. [ASSUMED] a guess, swept.
# [LITERATURE] from a paper. [PROPERTY] physical constant. [TARGET] a goal. [MEASURED] none yet.
import numpy as np
from dataclasses import dataclass

# water at ~22 C
RHO = 998.0          # kg/m^3   [PROPERTY]
MU = 9.55e-4         # Pa s     [PROPERTY]
NU = MU / RHO        # m^2/s    [PROPERTY]
G = 9.81             # m/s^2    [PROPERTY]


@dataclass(frozen=True)
class Jar:
    D: float = 0.0744           # m inner diameter       [DESIGN]
    H: float = 0.0737           # m liquid height        [DESIGN]
    Q_air_Lpm: float = 0.25     # L/min pump             [DESIGN]
    y_outlet: float = 0.058     # m A outlet holes       [DESIGN]
    x_stone_frac: float = 0.35  # B stone, fraction of W [DESIGN]
    V_liquid_mL: float = 320.0  # mL                     [DESIGN]

    @property
    def area(self):
        return np.pi * self.D ** 2 / 4

    @property
    def U_G(self):
        # superficial gas velocity
        return (self.Q_air_Lpm / 1000 / 60) / self.area


JAR = Jar()

# solver sweep
U0_SWEEP = (0.02, 0.05, 0.12)       # m/s     [ASSUMED]
DT_SWEEP = (1e-5, 3e-5, 1e-4)       # m^2/s   [ASSUMED]
STAGNANT_FACTOR = 0.1               # D_t cut above A outlets  [ASSUMED]

# biology time scales, seconds (other strains)
DAY, HOUR = 86400.0, 3600.0         # s       [PROPERTY]
BIO_TIMESCALES = {
    "bead release": (2 * DAY, 7 * DAY),
    "root capture": (1 * HOUR, 48 * HOUR),
    "die-off": (1.3 * DAY, 20 * DAY),
}
