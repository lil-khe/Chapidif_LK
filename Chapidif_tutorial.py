"""
Chapidif tutorial script
=========================
Nothing here is required beyond what's marked "REQUIRED".

Sections:
  1. Imports and dummy app
  2. Material
  3. Dielectric-function model choice
  4. Oscillators (Mermin/DL/Drude/...)
  5. GOS (core-shell ionization edges)
  6. Beam / calculation-quality parameters
  7. Exchange-correction settings
  8. MELF mode (per-oscillator collective vs. individual excitation, Egap)
  9. initParArray() -- always last
  10. Plots
"""

import numpy as np
from calculate import calculate
from run_and_plot import run_and_plot
from read_write_files import read_write_files
from null_app import NullApp                 # dummy app so Chapidif is happy outside its own GUI
from plots import *


# ---------------------------------------------------------------------------
# 1. Imports and dummy app
# ---------------------------------------------------------------------------
calc = calculate(parent=NullApp)


def addMermin(amp, omega, gamma, i, Eth=0.0, Delta=0.0):
    """Fill oscillator slot i."""
    calc.Amps[i] = amp        # REQUIRED. Oscillator strength/amplitude.
    calc.Omegas[i] = omega    # REQUIRED. Oscillator position (eV).
    calc.Gammas[i] = gamma    # REQUIRED. Damping/width (eV).
    calc.Eth[i] = Eth         # optional. Edge-factor threshold energy (eV): F(E-Eth)=1/(1+exp(-Delta*(E-Eth))).
    calc.Delta[i] = Delta     # optional. Edge-factor sharpness (eV^-1). Delta=0 -> factor is always 1 (feature off).


def addGOS(occ, edge, nl, Z, i):
    """Fill GOS (Generalized Oscillator Strength) slot i, for core-shell ionization edges."""
    calc.ConcGOS[i] = occ     # REQUIRED. Number of electrons in this shell per unit cell.
    calc.EdgeGOS[i] = edge    # REQUIRED. Binding/edge energy of the shell (eV).
    calc.nlGOS[i] = nl        # REQUIRED. Shell code: 10=1s, 20=2s, 21=2p, 30=3s, 31=3p, 32=3d.
    calc.ZGOS[i] = Z          # REQUIRED. Atomic number of the element this shell belongs to.


# ---------------------------------------------------------------------------
# 2. Material
# ---------------------------------------------------------------------------
calc.specificweight = 2.329     # REQUIRED. Mass density (g/cm^3).
calc.massunitcell = 28.085     # REQUIRED. Molar mass of the formula unit (g/mol). Together with specificweight this fixes calc.UnitCellDensity (per Angstrom^3).

# calc.InitDF() zeroes all oscillator arrays: call it once before you start filling anything in, and never again afterwards.
calc.InitDF()
calc.get_df_properties()      # optional here, but useful to sanity check densities/sum rules as you go


# ---------------------------------------------------------------------------
# 3. Dielectric-function model choice
# ---------------------------------------------------------------------------
# calc.DFmodel selects which physics the oscillators represent
calc.DFmodel = "Mermin"


# ---------------------------------------------------------------------------
# 4. Oscillators
# ---------------------------------------------------------------------------
# Slot index i must be unique per oscillator, 0 <= i < calc.maxOscillators (250).
# Eth/Delta are optional per-oscillator "switch-on" parameters: leave both at 0.0 for ordinary behaviour (F=1 everywhere, i.e. no edge suppression).

addMermin(0.96, 16.7, 3, i=0, Eth=1.2, Delta=1)			# main valence plasmon
addMermin(4.5e-2, 120, 120, i=1, Eth=100, Delta=1000000)	# inner-shell-adjacent oscillator; Delta->inf reproduces a step function


# ---------------------------------------------------------------------------
# 5. GOS (core-shell ionization edges)
# ---------------------------------------------------------------------------
calc.ApplySumRuleToGOS = 1    # 1: rescale GOS with q so the f-sum rule is respected at all q; 0: off
addGOS(occ=2, edge=1839, nl=10, Z=14, i=0)			# Si 1s; Energy comes from https://winter.group.shef.ac.uk/webelements/silicon/atoms.html


# ---------------------------------------------------------------------------
# 6. Beam / calculation-quality parameters
# ---------------------------------------------------------------------------
calc.particle = "electron"        # "electron" or "proton"
calc.E0 = 50.0                    # beam energy (keV): used by eps/ELF plots (q, energy-loss axis);
calc.q = 0                        # momentum transfer (a.u.) used for single-q plots (ELF, eps, dispersion slice)

calc.LowerELimit = 0              # eV, start of the energy-loss axis for eps/ELF plots
calc.UpperELimit = 100            # eV, end of the energy-loss axis for eps/ELF plots
calc.Stepsize = 0.1               # eV, step of that axis

calc.maxEnergyDensityEffect = 1e10 # eV; energy below which the density-effect correction to the real part is evaluated. Slow but needed for high-quality results.

calc.Stopping_calc_quality = 1    # 0 (fast), 1, or 2 (slow/most precise): controls the q-integration precision and the low-energy-tail momentum grid used in calc.calccurves() (IMFP/stopping/straggling/cross-section).

calc.MottCorrection = 1           # protons only: 1 applies the Mott correction to DIIMFP
calc.Dispersion_relativistic = 1  # 1: relativistic q-omega dispersion; 0: non-relativistic


# ---------------------------------------------------------------------------
# 7. Exchange-correction settings (electrons only)
# ---------------------------------------------------------------------------
#Pick one approach only.
calc.AddELF = 1                    # Required to be 1 for the Eth/Delta edge factor to act exactly on each oscillator's ELF.
calc.ExchangeCorrection = False    # False: no exchange correction (fastest); True: Ashley's method by default (or SBethe's, see Exchange_as_in_SBethe below)
calc.Exchange_as_in_SBethe = False # only used if ExchangeCorrection=True: False=Ashley, True=SBethe
calc.BE_for_exchange = 50.0        # eV, mean binding energy used by Ashley/SBethe exchange
calc.UseBornOchkurExchange = False # alternative, simpler, NON-relativistic exchange factor (1 - x + x^2, x=q^2/2E0), applied globally to the whole ELF instead of Ashley/SBethe.


# ---------------------------------------------------------------------------
# 8. MELF mode (optional): per-oscillator collective vs. individual excitation
# ---------------------------------------------------------------------------
# UseMELF=True switches calc.calccurves()/DIIMFP from "one combined eps" to "sum of independent per-oscillator ELF_i" (de Vera & Garcia-Molina, JPCC 2019, 123, 2075), letting each Mermin oscillator be either:
#   - a COLLECTIVE excitation (a plasmon): no exchange, integrated up to E0
#   - an INDIVIDUAL excitation (an ionization): Born-Ochkur exchange applied, integrated only up to (E0+BindingEnergy)/2
# GOS shells are always treated as individual/ionization automatically in this mode: nothing to set for them.
calc.UseMELF = False

# Egap: lower integration bound (eV): energy losses below the gap are skipped entirely.
calc.Egap = 1.12

# Per-oscillator switches, one entry per oscillator index i:
#   ExcitationTypeForOsc[i] = 1  -> collective (plasmon): UseExchangeForOsc should be False
#   ExcitationTypeForOsc[i] = 2  -> individual (ionization): UseExchangeForOsc should be True
calc.ExcitationTypeForOsc[0] = 1       # oscillator 0: collective
calc.UseExchangeForOsc[0] = False      # oscillator 0: no exchange
calc.ExcitationTypeForOsc[1] = 1       # oscillator 1: individual
calc.UseExchangeForOsc[1] = True       # oscillator 1: exchange on
calc.BindingEnergyForOsc[1] = 120      # oscillator 1: binding energy (required if exchange on)


# ---------------------------------------------------------------------------
# 9. initParArray(): must be called after ALL of the above, and again any time you change a value afterwards
# ---------------------------------------------------------------------------
calc.initParArray()


# ---------------------------------------------------------------------------
# 10. Plots
# ---------------------------------------------------------------------------
# All plotting helpers live in plots.py (imported above.
#
# General pattern: every function takes the 'calc' object, does whatever calculation is needed internally, and returns (fig, ax) so you can keep tweaking the plot afterwards (titles, limits, saving under a different name, etc.). Every function also accepts an optional custom 'label'/'color' (omit them to just use the built-in defaults). For some plots you can set loga=True.
#
# Figures are shown interactively and simultaneously saved as PNGs under ./FIT/ (can be modified in plots.py).

# --- 10.1 Single dielectric-function quantities at the current calc.q ---
# These use calc.q and the [LowerELimit, UpperELimit] energy-loss window
# set in section 6. Change calc.q before calling if you want a different
# momentum transfer.

#plot_eps(calc)                     # Re[eps(omega)] and Im[eps(omega)]
#plot_ELF(calc)                     # energy loss function Im[-1/eps(omega)]
#plot_kk_test(calc)                 # Kramers-Kronig self-consistency check (Re/Im[eps] vs. their KK reconstruction)


# --- 10.2 Sum rules ---
# Sanity checks: at high enough omega, the F-sum / Bethe-sum rules should
# converge towards the number of valence (+ core, if GOS included)
# electrons per unit cell; the KK sum rule should converge towards 1.

#plot_Z_sum_rules(calc)
#plot_kk_sum_rules(calc)


# --- 10.3 IMFP / stopping power / straggling / cross section vs. beam energy ---
# These call calc.calccurves() internally (recompute=True by default),
# which scans projectile energy according to calc.NStopping / calc.IncrFactor
# / calc.first_electron_energy (or first_proton_energy) unless you override the
# window with e_min_keV/e_max_keV below.

#plot_imfp(calc, overlay_tpp=True)   # IMFP vs. energy; overlay_tpp adds the TPP-2M estimate (electrons only)
#plot_stopping(calc)
#plot_straggling(calc)
#plot_cross_section(calc)

# Or all three stacked in one figure, sharing the x-axis:
#plot_all_three(calc)

# Restrict to a specific energy window instead of whatever NStopping/IncrFactor
# currently produce (this overwrites calc.NStopping/IncrFactor/first_*_energy):
#plot_imfp(calc, e_min_keV=1e-1, e_max_keV=5e1, n_points=80)


# --- 10.4 2D q-omega maps ---
# Plot the ELF dispersion with energy and momentum.

#plot_dispersion(calc, q_max=5.0, q_step=0.03)


# --- 10.5 Custom labels/colors and overlaying multiple curves on one axis ---
# Pass an existing 'ax' to draw several calc configurations on the same
# plot (e.g. comparing two materials, two oscillators/GOS settings or with
# experimental reference side by side).

#fig, ax = plt.subplots(figsize=(6, 6))
#data_exp = np.loadtxt("Reference.txt")
#energy = data_exp[:,0]
#intensity = data_exp[:,1]
#ax.scatter(energy, intensity, color="black", label="Reference")
#plot_ELF(calc, ax=ax, label="Si MELF", color="blue")
# ... change oscillators / re-run initParArray() here ...
# plot_ELF(calc, ax=ax, label="Si MELF + GOS", color="orange")
#ax.legend()
#plt.show()




