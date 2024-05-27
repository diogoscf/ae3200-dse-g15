import numpy as np
import json
import os
import matplotlib.pyplot as plt

import sys

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from isa import isa

FILE = "conventional.json"
COMMUTER = False

ft_to_m = 0.3048
FL_to_m = 0.3048 * 100.0
m_to_ft = 1 / 0.3048
m_to_FL = 1 / (0.3048 * 100.0)
lbs_to_kg = 0.45359237
kg_to_lbs = 1 / 0.45359237

G = 9.80665  # [m/s^2]

FONTSIZE = 40


if __name__ == "__main__":
    aircraft_data = json.load(open(os.path.join(os.path.dirname(__file__), "..", "Configurations", FILE), "r"))
    # MTOW_lbs = aircraft_data["MTOW_lbs"]
    # MTOW_kg = MTOW_lbs * lbs_to_kg
    MTOW_N = aircraft_data["MTOW_N"]
    MTOW_kg = MTOW_N / G

    if MTOW_kg > 8618:
        raise ValueError(
            "MTOW exceeds 8618 kg, so this is not a CS-23 aircraft. This module is only for CS-23 aircraft."
        )

    h = 0  # [m]
    t, p, rho = isa(h)
    t_0, p_0, rho_0 = isa(0)

    WS_kgm2 = aircraft_data["W/S_N/m2"]
    WS_lbft2 = aircraft_data["W/S_N/m2"] * kg_to_lbs / (m_to_ft**2)
    print(f"W/S: {WS_kgm2} kg/m^2 ({WS_lbft2:.3f} lb/ft^2)")

    MGC_m = aircraft_data["MGC_m"]
    clalpha = aircraft_data["CLalpha"]

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

    V_S0 = np.sqrt(2 * WS_kgm2 / (rho * CLmax_land))
    V_S1 = np.sqrt(2 * WS_kgm2 / (rho * CLmax_clean))

    mu_g = 2 * WS_kgm2 / (rho * MGC_m * clalpha * G)
    k_g = 0.88 * mu_g / (5.3 + mu_g)
    print(k_g)

    Ude_cruise_fps = 50
    Ude_cruise_ms = Ude_cruise_fps * 0.3048

    Ude_dive_fps = 25
    Ude_dive_ms = Ude_dive_fps * 0.3048

    Ude_B_fps = 66  # Only for commuter airplanes
    Ude_B_ms = Ude_B_fps * 0.3048

    n_cruise_pve = 1 + (k_g * rho_0 * Ude_cruise_ms * Vc_ms * clalpha) / (2 * WS_kgm2)
    n_cruise_nve = 1 - (k_g * rho_0 * Ude_cruise_ms * Vc_ms * clalpha) / (2 * WS_kgm2)

    n_dive_pve = 1 + (k_g * rho_0 * Ude_dive_ms * Vd_ms * clalpha) / (2 * WS_kgm2)
    n_dive_nve = 1 - (k_g * rho_0 * Ude_dive_ms * Vd_ms * clalpha) / (2 * WS_kgm2)

    V_B_intersect = (
        k_g * rho_0 * Ude_B_ms * clalpha
        + np.sqrt((k_g * rho_0 * Ude_B_ms * clalpha) ** 2 + 8 * WS_kgm2 * rho * CLmax_clean)
    ) / (2 * rho * CLmax_clean)
    V_B_stall = V_S1 * np.sqrt(n_cruise_pve)

    V_B = min(V_B_intersect, V_B_stall)

    n_B_pve = 1 + (k_g * rho_0 * Ude_B_ms * V_B * clalpha) / (2 * WS_kgm2)  # Commuter only
    n_B_nve = 1 - (k_g * rho_0 * Ude_B_ms * V_B * clalpha) / (2 * WS_kgm2)

    n_max = max(n_B_pve, n_cruise_pve, n_dive_pve) if COMMUTER else max(n_cruise_pve, n_dive_pve)
    n_min = min(n_B_nve, n_cruise_nve, n_dive_nve) if COMMUTER else min(n_cruise_nve, n_dive_nve)

    print(f"n_max: {n_max:.2f}")
    print(f"n_min: {n_min:.2f}")

    # Gust Diagram
    plt.rcParams.update({"font.size": FONTSIZE})
    fig, ax = plt.subplots()
    fig.tight_layout()
    ax.axhline(linewidth=2, color="k")
    ax.axvline(linewidth=2, color="k")
    style = "r-"
    LIGHTCOLOUR = "dimgrey"

    if COMMUTER:
        ax.plot([0, V_B], [1, n_B_pve], style, linewidth=2, zorder=20)
        ax.plot([0, Vc_ms], [1, n_cruise_pve], color=LIGHTCOLOUR, linestyle="--", linewidth=2, zorder=20)
        ax.plot([V_B, Vc_ms], [n_B_pve, n_cruise_pve], style, linewidth=2, zorder=20)

        ax.plot([0, V_B], [1, n_B_nve], style)
        ax.plot([0, Vc_ms], [1, n_cruise_nve], color=LIGHTCOLOUR, linestyle="--", linewidth=2, zorder=20)
        ax.plot([V_B, Vc_ms], [n_B_nve, n_cruise_nve], style, linewidth=2, zorder=20)
    else:
        ax.plot([0, Vc_ms], [1, n_cruise_pve], style, linewidth=2, zorder=20)
        ax.plot([0, Vc_ms], [1, n_cruise_nve], style, linewidth=2, zorder=20)

    ax.plot([0, Vd_ms], [1, n_dive_pve], color=LIGHTCOLOUR, linestyle="--", linewidth=2, zorder=20)
    ax.plot([0, Vd_ms], [1, n_dive_nve], color=LIGHTCOLOUR, linestyle="--", linewidth=2, zorder=20)

    ax.plot([Vc_ms, Vd_ms], [n_cruise_pve, n_dive_pve], style, linewidth=2, zorder=20)
    ax.plot([Vc_ms, Vd_ms], [n_cruise_nve, n_dive_nve], style, linewidth=2, zorder=20)

    ax.plot([Vd_ms, Vd_ms], [n_dive_nve, n_dive_pve], style, linewidth=2, zorder=20)

    markersize = 15
    if COMMUTER:
        ax.plot(V_B, n_B_pve, "ro", ms=markersize, zorder=20)
        # ax.annotate("B'", (V_B, n_B_pve), textcoords="offset points", xytext=(-2, 8), ha="center", fontweight="bold", fontsize=15)
        ax.plot(V_B, n_B_nve, "ro", ms=markersize, zorder=20)
        # ax.annotate("G'", (V_B, n_B_nve), textcoords="offset points", xytext=(-2, -20), ha="center", fontweight="bold", fontsize=15)

    ax.plot(Vc_ms, n_cruise_pve, "ro", ms=markersize, zorder=20)
    # ax.annotate("C'", (Vc_ms, n_cruise_pve), textcoords="offset points", xytext=(2, 8), ha="center", fontweight="bold", fontsize=15)
    ax.plot(Vc_ms, n_cruise_nve, "ro", ms=markersize, zorder=20)
    # ax.annotate("F'", (Vc_ms, n_cruise_nve), textcoords="offset points", xytext=(2, -20), ha="center", fontweight="bold", fontsize=15)

    ax.plot(Vd_ms, n_dive_pve, "ro", ms=markersize, zorder=20)
    # ax.annotate("D'", (Vd_ms, n_dive_pve), textcoords="offset points", xytext=(4, 8), ha="center", fontweight="bold", fontsize=15)
    ax.plot(Vd_ms, n_dive_nve, "ro", ms=markersize, zorder=20)
    # ax.annotate("E'", (Vd_ms, n_dive_nve), textcoords="offset points", xytext=(4, -20), ha="center", fontweight="bold", fontsize=15)

    ax.plot(0, 1, "ro", ms=markersize, zorder=20)

    if COMMUTER:
        ax.plot([V_B, V_B], [n_B_nve, n_B_pve], linestyle="--", linewidth=2, color=LIGHTCOLOUR, zorder=-1)
        # ax.plot(Vc_ms, 0, marker="o", color="grey")
        ax.annotate(
            "$V_B$",
            (V_B, 0),
            textcoords="offset points",
            xytext=(30, 10),
            ha="center",
            fontsize=FONTSIZE,
            color=LIGHTCOLOUR,
        )

    ax.plot([Vc_ms, Vc_ms], [n_cruise_nve, n_cruise_pve], linestyle="--", linewidth=2, color=LIGHTCOLOUR, zorder=-1)
    # ax.plot(Vc_ms, 0, marker="o", color="grey")
    ax.annotate(
        "$V_C$",
        (Vc_ms, 0),
        textcoords="offset points",
        xytext=(30, 10),
        ha="center",
        fontsize=FONTSIZE,
        color=LIGHTCOLOUR,
    )

    ax.annotate(
        "$V_D$",
        (Vd_ms, 0),
        textcoords="offset points",
        xytext=(-30, 10),
        ha="center",
        fontsize=FONTSIZE,
        color=LIGHTCOLOUR,
    )

    ax.plot([0, Vd_ms], [1, 1], linestyle="--", linewidth=2, color=LIGHTCOLOUR, zorder=-1)

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
        os.path.join(os.path.dirname(__file__), "..", "..", "Figures", f"vn-gust-{aircraft_data['name']}.pdf"),
        bbox_inches="tight",
        dpi=200,
    )
    plt.show()
