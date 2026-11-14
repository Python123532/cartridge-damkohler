from dataclasses import dataclass
import numpy as np

# ---- fluid (water, ~22 C) ---------------------------------------------------
RHO = 998.0          # kg/m^3
MU = 9.55e-4         # Pa s
NU = MU / RHO        # m^2/s
G = 9.81             # m/s^2


@dataclass(frozen=True)
class Jar:
    D: float = 0.0744          # m, inner diameter (CAD)
    H: float = 0.0737          # m, liquid height (CAD)
    Q_air_Lpm: float = 0.25    # L/min, pump setting
    y_outlet: float = 0.058    # m, Design A airlift outlet-hole height (CAD)
    x_stone_frac: float = 0.35  # Design B air stone offset, fraction of width
    V_liquid_mL: float = 320.0

    @property
    def area(self):
        return np.pi * self.D ** 2 / 4

    @property
    def U_G(self):
        """Superficial gas velocity, m/s."""
        return (self.Q_air_Lpm / 1000 / 60) / self.area


JAR = Jar()

# ---- biological time scales (seconds), ranges from the handoff -------------
DAY, HOUR = 86400.0, 3600.0
BIO_TIMESCALES = {
    "bead release": (2 * DAY, 7 * DAY),
    "root capture": (1 * HOUR, 48 * HOUR),
    "die-off": (1.3 * DAY, 20 * DAY),
}

# ---- bead / release ---------------------------------------------------------
BEAD_D_RANGE = (2e-3, 4e-3)       # m (Balla 2025 fresh beads 3.17 +/- 0.18 mm)
M_EXPONENT_RANGE = (0.4, 0.8)     # brackets Rittmann 1982's 0.58

# Local slip velocity past a bead, as a fraction of the jar circulation velocity
# V_c. ASSUMPTION: Design A beads sit in the
# airlift riser and see a large fraction of the riser liquid velocity; Design B
# beads sit in a sheltered basket and see a small fraction of ambient flow.
# These priors reproduce the handoff's tau ranges. They are the weakest inputs
# to the shear model and the obvious thing to measure (PIV / particle tracking).
SLIP_FRAC_A = (0.10, 0.50)
SLIP_FRAC_B = (0.005, 0.05)

# ---- dose / pre-registration -----------------------------------------------
# Balla et al. 2025, AIMS Microbiol 11(1):87-125, doi 10.3934/microbiol.2025006
BALLA_B25_PLATEAU_CFU_per_mL = 8.1e6   # day 60, 1 g beads in 10 mL
BALLA_TUBE_VOLUME_mL = 10.0
BALLA_PA_D1 = 1.38e7                   # Pantoea Pa, CFU/mL, day 1
BALLA_PA_D25 = 7.0e7                   # day 25
TARGET_BAND = (15_600, 31_200)         # CFU/mL
TARGET_DESIGN = 20_000
OVERSHOOT_CAP = 62_400
DIEOFF_K = (0.0, 0.05, 0.3, 0.8)       # 1/day
