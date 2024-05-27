import numpy as np
import json
import os
import matplotlib.pyplot as plt

import sys

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from isa import isa

FILE = "conventional.json"
PLOT_FLAPPED = False

FONT_SIZE = 40


def nmax_f(MTOW_lbs):
    return min(2.1 + 24000 / (MTOW_lbs + 10000), 3.8)


def nmin_f(MTOW_lbs):
    return -0.4 * nmax_f(MTOW_lbs)


def nult_f(MTOW_lbs):
    return 1.5 * nmax_f(MTOW_lbs)


ft_to_m = 0.3048
FL_to_m = 0.3048 * 100.0
m_to_ft = 1 / 0.3048
m_to_FL = 1 / (0.3048 * 100.0)
lbs_to_kg = 0.45359237
kg_to_lbs = 1 / 0.45359237

G = 9.80665  # [m/s^2]


if __name__ == "__main__":
    aircraft_data = json.load(open(os.path.join(os.path.dirname(__file__), "..", "Configurations", FILE), "r"))
    # MTOW_lbs = aircraft_data["MTOW_lbs"]
    # MTOW_kg = MTOW_lbs * lbs_to_kg
    MTOW_N = aircraft_data["MTOW_N"]
    MTOW_kg = MTOW_N / G
    MTOW_lbs = MTOW_kg * kg_to_lbs
    print(MTOW_lbs)

    if MTOW_kg > 8618:
        raise ValueError(
            "MTOW exceeds 8618 kg, so this is not a CS-23 aircraft. This module is only for CS-23 aircraft."
        )

    nmax = nmax_f(MTOW_lbs)
    nmin = nmin_f(MTOW_lbs)

    h = 0  # [m]
    t, p, rho = isa(h)

    WS_kgm2 = aircraft_data["W/S_N/m2"]
    WS_lbft2 = aircraft_data["W/S_N/m2"] * kg_to_lbs / (m_to_ft**2)
    print(f"W/S: {WS_kgm2} kg/m^2 ({WS_lbft2:.3f} lb/ft^2)")

    Vc_ms = aircraft_data["Vc_m/s"]

    # CS-23.335(b)(2)(i)
    k_vc_vd = 1.4
    if WS_lbft2 > 20:
        k_vc_vd = np.interp(WS_lbft2, [20, 100], [1.4, 1.35])
    Vd_ms = k_vc_vd * Vc_ms
    print(f"V_D: {Vd_ms:.2f} m/s")

    V_ms = np.linspace(0, np.ceil(Vd_ms), 1000)
    q = 0.5 * rho * V_ms**2

    CLmax_clean = aircraft_data["CLmax_clean"]
    CLmax_land = aircraft_data["CLmax_land"]

    V_A = np.sqrt(2 * nmax * WS_kgm2 / (rho * CLmax_clean))
    V_S1 = np.sqrt(2 * WS_kgm2 / (rho * CLmax_clean))
    V_H = np.sqrt(2 * (-nmin) * WS_kgm2 / (rho * CLmax_clean))

    if V_A > Vc_ms:
        raise ValueError("V_A > V_C")

    V_S0 = np.sqrt(2 * WS_kgm2 / (rho * CLmax_land))

    n_OA = q * CLmax_clean / WS_kgm2
    n_OA = np.concatenate((n_OA[np.where(n_OA < nmax)], [nmax]))
    V_OA = np.concatenate((V_ms[: np.shape(n_OA)[0] - 1], [V_A]))

    nmax_flap = 2  # TODO: Check
    V_I = np.sqrt(2 * 2 * WS_kgm2 / (rho * CLmax_land))
    n_OI = q * CLmax_land / WS_kgm2
    n_OI = np.concatenate((n_OI[np.where(n_OI < nmax_flap)], [nmax_flap]))  # needs to be changed
    V_OI = np.concatenate((V_ms[: np.shape(n_OI)[0] - 1], [V_I]))

    V_J = np.sqrt(2 * nmax_flap * WS_kgm2 / (rho * CLmax_clean))
    n_IJ = [nmax_flap, nmax_flap]
    V_IJ = [V_I, V_J]

    n_AD = [nmax, nmax]
    V_AD = [V_A, Vd_ms]

    n_OH = -n_OA[np.where(n_OA < -nmin)]
    n_OH = np.concatenate((n_OH, [nmin]))
    V_OH = np.concatenate((V_ms[: np.shape(n_OH)[0] - 1], [V_H]))

    n_HF = [nmin, nmin]
    V_HF = [V_H, Vc_ms]

    n_FE = [nmin, 0]
    V_FE = [Vc_ms, Vd_ms]

    n_ED = [0, nmax]
    V_ED = [Vd_ms, Vd_ms]

    print("V_S1:", V_S1)
    print("V_S0:", V_S0)
    print("nmax:", nmax)
    print("nmin:", nmin)

    plt.rcParams.update({"font.size": FONT_SIZE})
    fig, ax = plt.subplots()
    fig.tight_layout()
    ax.axhline(linewidth=2, color="k")
    ax.axvline(linewidth=2, color="k")

    style = "r-"
    ax.plot(V_OA, n_OA, style, linewidth=2, zorder=20)
    ax.plot(V_AD, n_AD, style, linewidth=2, zorder=20)
    ax.plot(V_OH, n_OH, style, linewidth=2, zorder=20)
    ax.plot(V_HF, n_HF, style, linewidth=2, zorder=20)
    ax.plot(V_FE, n_FE, style, linewidth=2, zorder=20)
    ax.plot(V_ED, n_ED, style, linewidth=2, zorder=20)

    markersize = 15
    if PLOT_FLAPPED:
        ax.plot(V_OI, n_OI, style, linewidth=2, zorder=20)
        ax.plot(V_IJ, n_IJ, style, linewidth=2, zorder=20)

    ax.plot(V_A, nmax, "ro", ms=markersize)
    # ax.annotate("A", (V_A, nmax), textcoords="offset points", xytext=(-2, 8), ha="center", fontweight="bold", fontsize=FONT_SIZE)

    ax.plot(Vd_ms, nmax, "ro", ms=markersize)
    # ax.annotate("D", (Vd_ms, nmax), textcoords="offset points", xytext=(2, 8), ha="center", fontweight="bold", fontsize=FONT_SIZE)

    ax.plot(Vd_ms, 0, "ro", ms=markersize)
    # ax.annotate("E", (Vd_ms, 0), textcoords="offset points", xytext=(8, 8), ha="center", fontweight="bold", fontsize=FONT_SIZE)

    ax.plot(Vc_ms, nmin, "ro", ms=markersize)
    # ax.annotate("F", (Vc_ms, nmin), textcoords="offset points", xytext=(2, -20), ha="center", fontweight="bold", fontsize=FONT_SIZE)

    ax.plot(V_H, nmin, "ro", ms=markersize)
    # ax.annotate("H", (V_H, nmin), textcoords="offset points", xytext=(-2, -20), ha="center", fontweight="bold", fontsize=FONT_SIZE)

    ax.plot(0, 0, "ro", ms=markersize)

    if PLOT_FLAPPED:
        ax.plot(V_I, nmax_flap, "ro", ms=markersize)
        # ax.annotate("I", (V_I, nmax_flap), textcoords="offset points", xytext=(-2, 8), ha="center", fontweight="bold", fontsize=15)

        ax.plot(V_J, nmax_flap, "ro", ms=markersize)
        # ax.annotate("J", (V_J, nmax_flap), textcoords="offset points", xytext=(-2, 8), ha="center", fontweight="bold", fontsize=15)

    LIGHTCOLOUR = "dimgrey"

    ax.plot([Vc_ms, Vc_ms], [nmin, 0], linewidth=2, linestyle="--", color=LIGHTCOLOUR, zorder=-1)
    # ax.plot(Vc_ms, 0, marker="o", color="grey")
    ax.annotate(
        "$V_C$",
        (Vc_ms, 0),
        textcoords="offset points",
        xytext=(0, 11),
        ha="center",
        fontsize=FONT_SIZE,
        color=LIGHTCOLOUR,
    )

    ax.annotate(
        "$V_D$",
        (Vd_ms, 0),
        textcoords="offset points",
        xytext=(-38, 11),
        ha="center",
        fontsize=FONT_SIZE,
        color=LIGHTCOLOUR,
    )

    ax.plot([V_A, V_A], [0, nmax], linewidth=2, linestyle="--", color=LIGHTCOLOUR, zorder=-1)
    # ax.plot(Vc_ms, 0, marker="o", color="grey")
    ax.annotate(
        "$V_A$",
        (V_A, 0),
        textcoords="offset points",
        xytext=(0, -38),
        ha="center",
        fontsize=FONT_SIZE,
        color=LIGHTCOLOUR,
    )

    ax.plot([V_S1, V_S1], [-1, 1], linewidth=2, linestyle="--", color=LIGHTCOLOUR, zorder=-1)
    # ax.plot(V_S1, 0, marker="o", color="grey")
    ax.annotate(
        "$V_{S1}$",
        (V_S1, 0),
        textcoords="offset points",
        xytext=(40, 11),
        ha="center",
        fontsize=FONT_SIZE,
        color=LIGHTCOLOUR,
    )

    ax.annotate(
        f"$n_{{max}} = {nmax:.2f}$",
        (V_A + (Vd_ms - V_A) / 2, nmax),
        textcoords="offset points",
        xytext=(0, -38),
        ha="center",
        fontsize=FONT_SIZE,
    )
    ax.annotate(
        f"$n_{{min}} = {nmin:.2f}$",
        (V_H + (Vc_ms - V_H) / 2, nmin),
        textcoords="offset points",
        xytext=(0, 10),
        ha="center",
        fontsize=FONT_SIZE,
    )

    ax.plot([0, Vd_ms], [1, 1], linewidth=2, linestyle="--", color=LIGHTCOLOUR, zorder=5)

    ax.set_xlabel("Velocity (EAS) [m/s]")
    ax.set_ylabel("Load Factor [-]")

    ax.set_ylim(-2.2, 4.2)
    ax.set_yticks(np.arange(-2.0, 4.2, 1))

    ax.set_xticks(np.arange(0, 90.1, 20))
    ax.set_xticks(np.arange(0, 90.1, 10), minor=True)

    ax.tick_params("both", length=10, width=1, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")

    ax.grid()

    fig.set_size_inches(16, 9)
    fig.tight_layout()
    fig.savefig(
        os.path.join(os.path.dirname(__file__), "..", "..", "Figures", f"vn-loading-{aircraft_data['name']}.pdf"),
        bbox_inches="tight",
        dpi=200,
    )
    plt.show()
