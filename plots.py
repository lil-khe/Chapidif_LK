import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Use an interactive backend so figures pop up and stay interactive
# instead of appearing only once saved to disk. Comment out if you
# run this in a notebook that already sets its own backend (%matplotlib widget/qt).
for _backend in ("TkAgg", "Qt5Agg", "QtAgg", "GTK3Agg", "WXAgg"):
    try:
        matplotlib.use(_backend, force=True)
        break
    except ImportError:
        continue
else:
    print("Warning: no interactive backend found, falling back to Agg (plots will only save, not display)")
plt.ion()


PATH = "./FIT/"
os.makedirs(PATH, exist_ok=True)

###	PLOT IMFP/STOPPING/STRAGGLING	###

_YLABELS = {
    "imfp": "IMFP (Å)",
    "stopping": "stopping (eV/Å)",
    "straggling": "straggling (eV²/Å)",
    "cross_section": "cross section (Å²)"
}
_ARRAY_NAMES = {
    "imfp": "IMFPEnergy",
    "stopping": "StoppingEnergy",
    "straggling": "StragglingEnergy",
    "cross_section": "CrosssectionEnergy"
}


def _save_and_show(fig, name):
    """Save a figure to PATH as a PNG and display it."""
    fig.savefig(os.path.join(PATH, name + ".png"), dpi=300, bbox_inches="tight")
    plt.show(block=True)


def set_energy_window(calc, e_min_keV, e_max_keV, n_points=50):
    """
    Set the projectile-energy window over which calc.calccurves() computes
    IMFP/stopping/straggling, by target min/max energy instead of raw
    IncrFactor. Mirrors the formula Chapidif itself uses internally
    (see calculate.py's Energy_Depth_Dist): IncrFactor is chosen so that
    NStopping steps of CurrentE *= IncrFactor go from e_min to e_max.
    """
    if n_points < 2:
        raise ValueError("n_points must be >= 2 to define an energy window")
    calc.NStopping = n_points
    calc.IncrFactor = (e_max_keV / e_min_keV) ** (1.0 / (n_points - 1))
    if calc.particle == "electron":
        calc.first_electron_energy = e_min_keV
    else:
        calc.first_proton_energy = e_min_keV


def plot_quantity(calc, quantity="imfp", x_axis_keV=True, recompute=True, loga=False,
                   overlay_tpp=False, ax=None, e_min_keV=None, e_max_keV=None,
                   n_points=None, label=None, color=None):
    """
    Plot one of "imfp" / "stopping" / "straggling" / "cross_section" vs.
    projectile energy (or velocity). This is the generic worker function
    used by plot_imfp / plot_stopping / plot_straggling / plot_cross_section.

    quantity      : which curve to plot, one of _ARRAY_NAMES keys
    x_axis_keV    : True -> x-axis is projectile energy (keV)
                    False -> x-axis is projectile velocity (a.u.)
    recompute     : call calc.calccurves(False) first (set False to reuse
                    curves you just computed elsewhere)
    loga          : log-log axes instead of linear
    overlay_tpp   : overlay the TPP-2M IMFP estimate (only meaningful for
                    quantity="imfp", electrons)
    ax            : plot into an existing axes instead of creating a new figure
    e_min_keV/e_max_keV/n_points : optionally set the energy window via
                    set_energy_window() before computing
    label, color  : optional custom legend label / line color; if omitted,
                    falls back to the default label and matplotlib's
                    automatic color cycling

    Returns (fig, ax) -- fig is None if `ax` was supplied.
    """
    if quantity not in _ARRAY_NAMES:
        raise ValueError(f"quantity must be one of {list(_ARRAY_NAMES)}, got {quantity!r}")

    if e_min_keV is not None and e_max_keV is not None:
        set_energy_window(calc, e_min_keV, e_max_keV, n_points or calc.NStopping)
        recompute = True
    elif (e_min_keV is None) != (e_max_keV is None):
        raise ValueError("e_min_keV and e_max_keV must be given together")

    if recompute:
        calc.calccurves(False)

    if x_axis_keV:
        x = calc.CurvesEnergy
        xlabel = calc.particle + " energy (keV)"
    else:
        x = calc.CurvesVelocity
        xlabel = calc.particle + " velocity (a.u.)"

    y = getattr(calc, _ARRAY_NAMES[quantity])
    plot_label = label if label is not None else _YLABELS[quantity]
    plot_kwargs = {"label": plot_label}
    if color is not None:
        plot_kwargs["color"] = color

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    ax.plot(x, y, **plot_kwargs)
    if loga:
        ax.set_xscale('log')
        ax.set_yscale('log')
        # bottom=0 is meaningless on a log axis, so only clip on linear plots
    else:
        ax.set_ylim(0, np.max(y) * 1.1)

    if overlay_tpp and quantity == "imfp" and calc.particle == "electron" and np.amax(calc.TPP_IMFPEnergy) > 0.0:
        ax.plot(x, calc.TPP_IMFPEnergy, color="purple", linestyle="dotted",
                label=f"TPP-2M ω_p={calc.w_p_TPP:.2f} eV, ρ={calc.specificweight:.2f} g/cm³")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(_YLABELS[quantity])
    ax.set_xlim(x[0], x[-1])
    ax.legend()

    if fig is not None:
        _save_and_show(fig, quantity)
    return fig, ax


def plot_imfp(calc, x_axis_keV=True, recompute=True, loga=False, overlay_tpp=False,
              e_min_keV=None, e_max_keV=None, n_points=None, label=None, color=None, ax=None):
    """Convenience wrapper: plot the inelastic mean free path (IMFP)."""
    return plot_quantity(calc, "imfp", x_axis_keV, recompute, loga, overlay_tpp,
                          e_min_keV=e_min_keV, e_max_keV=e_max_keV, n_points=n_points,
                          label=label, color=color, ax=ax)


def plot_stopping(calc, x_axis_keV=True, recompute=True, loga=False,
                   e_min_keV=None, e_max_keV=None, n_points=None, label=None, color=None, ax=None):
    """Convenience wrapper: plot the stopping power."""
    return plot_quantity(calc, "stopping", x_axis_keV, recompute, loga,
                          e_min_keV=e_min_keV, e_max_keV=e_max_keV, n_points=n_points,
                          label=label, color=color, ax=ax)


def plot_straggling(calc, x_axis_keV=True, recompute=True, loga=False,
                     e_min_keV=None, e_max_keV=None, n_points=None, label=None, color=None, ax=None):
    """Convenience wrapper: plot the energy-loss straggling."""
    return plot_quantity(calc, "straggling", x_axis_keV, recompute, loga,
                          e_min_keV=e_min_keV, e_max_keV=e_max_keV, n_points=n_points,
                          label=label, color=color, ax=ax)


def plot_cross_section(calc, x_axis_keV=True, recompute=True, loga=False,
                        e_min_keV=None, e_max_keV=None, n_points=None, label=None, color=None, ax=None):
    """Convenience wrapper: plot the total inelastic cross section."""
    return plot_quantity(calc, "cross_section", x_axis_keV, recompute, loga,
                          e_min_keV=e_min_keV, e_max_keV=e_max_keV, n_points=n_points,
                          label=label, color=color, ax=ax)


def plot_all_three(calc, x_axis_keV=True, recompute=True, loga=False,  e_min_keV=None, e_max_keV=None, n_points=None, labels=None, colors=None):
    """
    Three stacked subplots (IMFP / stopping / straggling), sharing the x-axis.

    labels, colors : optional lists of 3 entries, one per subplot
                      (order: imfp, stopping, straggling). Leave as None
                      to use the default labels/colors for all three.
    """
    if recompute:
        calc.calccurves(False)

    labels = labels or [None, None, None]
    colors = colors or [None, None, None]
    if len(labels) != 3 or len(colors) != 3:
        raise ValueError("labels and colors must each have exactly 3 entries")

    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(6, 10))
    for ax, quantity, lab, col in zip(axs, ("imfp", "stopping", "straggling"), labels, colors):
        plot_quantity(calc, quantity, x_axis_keV, recompute=False, loga=loga,
                      e_min_keV=e_min_keV, e_max_keV=e_max_keV, n_points=None, label=lab, color=col, ax=ax)
    axs[-1].set_xlabel(axs[0].get_xlabel())
    fig.subplots_adjust(hspace=0)
    _save_and_show(fig, "imfp_stopping_straggling")
    return fig, axs


###	PLOT FUNCTION	###

def plot_dispersion(calc, q_max=None, q_step=None, loga=False, max_eq=0.0, ax=None, fig=None, cmap="plasma"):
    """
    2D map of the loss function Im[-1/eps] vs energy loss omega (eV, y-axis)
    and momentum transfer q (a.u., x-axis), using Chapidif's own
    colorplot_lossfunction() + scale_image() (the same path run_and_plot's
    eq_plot() uses internally).

 
    q_max, q_step : override calc.UpperqLimit / calc.Stepsize_qplot before
                    plotting (omit to just use whatever calc is currently
                    set to). The q-axis always starts at 0.5*Stepsize_qplot,
                    matching Chapidif's own convention: there's no
                    separate q_min.
    loga          : log10 color scale (uses calc.LogXY / calc.scale_image's
                    built-in log handling) instead of linear.
    max_eq        : clip the color scale at this value (0.0 = auto-scale to
                    the data max, matching calc.max_eq's default meaning).
    """
    if q_max is not None:
        calc.UpperqLimit = q_max
    if q_step is not None:
        calc.Stepsize_qplot = q_step
    calc.LogXY = loga
    calc.max_eq = max_eq
    calc.bulk_eq = True  # bulk (not surface) energy-loss function, matches run_and_plot.eq_plot()
 
    # Changing UpperqLimit/Stepsize_qplot changes Nqstep/CenterFirstBin/etc,
    # which only get recomputed on initParArray() -- required any time a
    # setting changes after the initial setup (see tutorial section 9).
    calc.initParArray()
 
    calc.colorplot_lossfunction()
    calc.scale_image()
 
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
 
    ratio = calc.UpperqLimit / (calc.UpperELimit - calc.CenterFirstBin)
    im = ax.imshow(
        calc.my_scaled_image,
        extent=[0, calc.UpperqLimit, calc.UpperELimit, calc.LowerELimit],
        aspect=ratio,
        cmap=cmap
    )
    plt.colorbar(im, ax=ax, label="log10(loss function)" if loga else "loss function")

    ax.set_xlim(0, calc.UpperqLimit)
    ax.set_ylim(calc.CenterFirstBin, calc.UpperELimit)
    ax.set_xlabel("q (a.u.)")
    ax.set_ylabel("$\\omega$ (eV)")
 
    if fig is not None:
        _save_and_show(fig, "momentum_dispersion")
    return fig, ax


def plot_ELF(calc, loga=False, ax=None, label="ELF", color=None):
    """Plot the energy loss function Im[-1/eps(ω)] at the currently set calc.q."""
    calc.oneovereps1eps2()

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    plot_kwargs = {"label": label}
    if color is not None:
        plot_kwargs["color"] = color

    if loga:
        ax.loglog(calc.x_axis, calc.Result2, **plot_kwargs)
    else:
        ax.plot(calc.x_axis, calc.Result2, **plot_kwargs)

    ax.set_xlabel("$\\omega$ (eV)")
    ax.set_ylabel("Intensity (a.u.)")

    ax.set_xlim(calc.x_axis[0], calc.x_axis[-1])
    if loga:
        ax.set_ylim(1e-9, np.max(calc.Result2) * 10)
    else:
        ax.set_ylim(0, np.max(calc.Result2) * 1.1)

    if fig is not None:
        _save_and_show(fig, "ELF")
    return fig, ax


def plot_eps(calc, labels=None, colors=None, ax=None):
    """
    Plot Re[eps(ω)] and Im[eps(ω)] at the currently set calc.q.

    labels, colors : optional lists of 2 entries [Re label, Im label]
                      (and matching colors). Leave as None for defaults.
    """
    calc.eps1eps2()

    labels = labels or ["Re($\\varepsilon$)", "Im($\\varepsilon$)"]
    colors = colors or [None, None]
    if len(labels) != 2 or len(colors) != 2:
        raise ValueError("labels and colors must each have exactly 2 entries")

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    ax.plot(calc.x_axis, calc.Result1, label=labels[0], color=colors[0])
    ax.plot(calc.x_axis, calc.Result2, label=labels[1], color=colors[1])

    ax.set_xlabel("$\\omega$ (eV)")
    ax.set_ylabel("Intensity (a.u.)")

    ax.set_xlim(calc.x_axis[0], calc.x_axis[-1])
    ax.set_ylim(np.min(calc.Result1) * 1.1, max(np.max(calc.Result2), np.max(calc.Result1)) * 1.1)
    ax.legend()

    if fig is not None:
        _save_and_show(fig, "DF")
    return fig, ax


def plot_kk_test(calc, labels=None, colors=None, ax=None):
    """
    Kramers-Kronig consistency check: plot the computed Re[eps]/Im[eps]
    alongside their KK-transform reconstructions (Result3 = KK-transform
    of Im[eps], Result4 = KK-transform of Re[eps]). The two pairs should
    overlay closely if the dielectric function is KK-consistent.

    Returns (None, None) and prints a message if calc.eps_kk_test() fails.

    labels, colors : optional lists of 4 entries
                      [Re eps, Im eps, KK Im eps, KK Re eps]
    """
    error = calc.eps_kk_test()
    if error != 0:
        print("eps_kk_test failed, error code:", error)
        return None, None

    default_labels = ["Re($\\varepsilon$)", "Im($\\varepsilon$)",
                       "KK-transform Im($\\varepsilon$)", "KK-transform Re($\\varepsilon$)"]
    labels = labels or default_labels
    colors = colors or [None, None, None, None]
    if len(labels) != 4 or len(colors) != 4:
        raise ValueError("labels and colors must each have exactly 4 entries")

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    ax.plot(calc.x_axis, calc.Result1, label=labels[0], color=colors[0])
    ax.plot(calc.x_axis, calc.Result2, label=labels[1], color=colors[1])
    ax.plot(calc.x_axis, calc.Result3, linestyle="--", label=labels[2], color=colors[2])
    ax.plot(calc.x_axis, calc.Result4, linestyle="-.", label=labels[3], color=colors[3])

    ax.set_xlabel("$\\omega$ (eV)")
    ax.set_ylabel("Intensity (a.u.)")

    ax.set_xlim(calc.x_axis[0], calc.x_axis[-1])
    ax.legend()

    if fig is not None:
        _save_and_show(fig, "KK_sum_rule")
    return fig, ax


def plot_Z_sum_rules(calc, loga=False, ax=None, labels=None, colors=None):
    """
    Plot the F-sum rule, Bethe (-1/eps) sum rule, and k(ω) sum rule vs ω,
    each converging towards the number of electrons per unit cell.

    labels, colors : optional lists of 3 entries, one per curve
                      [F-sum, Bethe-sum, k-sum]
    """
    calc.sum_rules()

    default_labels = [
        r"$\frac{2}{\pi \Omega_p^2}\int_0^\omega \omega \, {\rm Im} [\epsilon (\omega, q=%s)] d\omega$ (F)" % str(calc.q),
        r"$\frac{2}{\pi \Omega_p^2}\int_0^\omega \omega \, {\rm Im} [-1/\epsilon (\omega, q=%s)] d\omega$ (Bethe)" % str(calc.q),
        r"$\frac{4}{\pi \Omega_p^2}\int_0^\omega \omega \, k(\omega,q=%s) d\omega$" % str(calc.q),
    ]
    labels = labels or default_labels
    colors = colors or [None, None, None]
    if len(labels) != 3 or len(colors) != 3:
        raise ValueError("labels and colors must each have exactly 3 entries")

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    ax.plot(calc.x_axis, calc.Result1, label=labels[0], color=colors[0])
    ax.plot(calc.x_axis, calc.Result2, label=labels[1], color=colors[1])
    ax.plot(calc.x_axis, calc.Result3, label=labels[2], color=colors[2])

    if loga:
        ax.set_xscale('log')

    ax.grid(color="grey")

    ax.set_xlabel("$\\omega$ (eV)")
    ax.set_ylabel("electrons per unit cell")

    ax.set_xlim(calc.x_axis[0], calc.x_axis[-1])
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
    ax.legend()

    ax.text(0.8, 0.85, "F at " + f"{calc.x_axis[-1]:.0f}" + " = " + f"{calc.Result1[-1]:.3f}",
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.8, 0.8, "Bethe at " + f"{calc.x_axis[-1]:.0f}" + " = " + f"{calc.Result2[-1]:.3f}",
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    if fig is not None:
        _save_and_show(fig, "Z_sum_rules")
    return fig, ax


def plot_kk_sum_rules(calc, loga=False, ax=None, label=None, color=None):
    """
    Plot the Kramers-Kronig (perfect-screening) sum rule, Result4, vs ω.

    label, color : optional custom legend label / line color
    """
    calc.sum_rules()

    default_label = r"$\frac{2}{\pi }\int_0^\omega \frac{1}{\omega}\, {\rm Im} [-1/\epsilon (\omega, q=%s)] d\omega$" % str(calc.q)
    plot_label = label if label is not None else default_label

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    ax.plot(calc.x_axis, calc.Result4, label=plot_label, color=color)

    if loga:
        ax.set_xscale('log')

    ax.set_xlabel("$\\omega$ (eV)")
    ax.set_ylabel("KK sum rule")

    ax.grid(color="grey")
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))

    ax.set_xlim(calc.x_axis[0], calc.x_axis[-1])
    ax.set_ylim(bottom=0)
    ax.legend()

    ax.text(0.8, 0.9, "= " + f"{calc.Result4[-1]:.3f}" + " at " + f"{calc.x_axis[-1]:.0f}",
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    if fig is not None:
        _save_and_show(fig, "amp_sum_rules")
    return fig, ax

###########################
