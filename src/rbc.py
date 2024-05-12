#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quick 'pocket calculator' for Rayleigh-Bénard Convection parameters. Comes
with a minimal command-line interface for quick inspection.

Provides classes ``RBC_cell()`` and ``RBC()`` useful for working out dataseries
on Rayleigh-Bénard Convection in your own scripts.
"""
__author__ = "Dennis van Gils"
__authoremail__ = "vangils.dennis@gmail.com"
__url__ = "https://github.com/Dennis-van-Gils/python-fluidprop"
__date__ = "13-05-2024"
__version__ = "1.2.0"

from typing import Union

import numpy as np
import numpy.typing as npt
import fluidprop
from fluidprop import pprint

# ------------------------------------------------------------------------------
#   RBC_Cell
# ------------------------------------------------------------------------------


class RBC_Cell:
    """Rayleigh-Bénard Convection cell.

    Args:
        name (`str`): Name of the cell.

        D (`float`): Cell diameter [m].

        L (`float`): Cell height [m].

        g (`float`): Local gravity [m/s^2]. Default 9.8 m/s^2.

    Properties:
        name (`str`): Name of the cell.

        D (`float`): Cell diameter [m].

        L (`float`): Cell height [m].

        g (`float`): Local gravity [m/s^2].

        area (`float`): Cell bottom plate area [m^2].

        Gamma (`float`): Aspect ratio of the cell: D / L.
    """

    def __init__(self, name: str, D: float, L: float, g: float = 9.8):
        self.name = name
        """Name of the cell"""
        self.D = D
        """Cell diameter [m]"""
        self.L = L
        """Cell height [m]"""
        self.g = g
        """Local gravity [m/s^2]"""
        self.area = np.pi * (D / 2) ** 2
        """Cell bottom plate area [m^2]"""
        self.Gamma = D / L
        """Aspect ratio of the cell: D / L"""


RBC_cells: list[RBC_Cell] = []
"""Predefined list of Rayleigh-Bénard Convection cells."""

# fmt: off
RBC_cells.append(RBC_Cell(name="unit cell 0.1 m"    , D=0.1  , L=0.1  , g=9.8))
RBC_cells.append(RBC_Cell(name="unit cell 1.0 m"    , D=1    , L=1    , g=9.8))
RBC_cells.append(RBC_Cell(name="HPCF-II, G=0.5"     , D=1.122, L=2.240, g=9.81))
RBC_cells.append(RBC_Cell(name="HPCF-IV, G=1"       , D=1.122, L=1.120, g=9.81))
RBC_cells.append(RBC_Cell(name="Barrel of Ilmenau"  , D=7.15 , L=6.3  , g=9.81))
RBC_cells.append(RBC_Cell(name="Chicago 1991, G=0.5", D=0.2  , L=0.4  , g=9.8))
RBC_cells.append(RBC_Cell(name="Chicago 1991, G=1"  , D=0.087, L=0.087, g=9.8))
RBC_cells.append(RBC_Cell(name="Chicago 1991, G=6.7", D=0.2  , L=0.03 , g=9.8))
# fmt: on

# ------------------------------------------------------------------------------
#   RBC
# ------------------------------------------------------------------------------


class RBC:
    """Container for Rayleigh-Bénard Convection parameters. Calculates
    non-dimensional parameters like the Rayleigh and Ekman number based on the
    given fluid, cell and the driving parameters `DT` and `Omega`.

    Args:
        fluid (`fluidprop.FluidProperties`):
            FluidProperties object evaluated at specific temperature(s) and
            pressure(s).

        cell (`rbc.RBC_Cell`):
            Rayleigh-Bénard Convection cell object.

        DT (`float` | `list[float]` | `numpy.ndarray[]`):
            Temperature difference between the plates [K].

        Omega (`float` | `list[float]` | `numpy.ndarray[]`, optional):
            Rotation rate along the cell axis [rad/s]. Default: 0 rad/s.

    """

    def __init__(
        self,
        fluid: fluidprop.FluidProperties,
        cell: RBC_Cell,
        DT: Union[float, list[float], npt.NDArray],
        Omega: Union[float, list[float], npt.NDArray] = 0,
    ):
        self.fluid = fluid
        self.cell = cell

        if isinstance(DT, (float, int)):
            DT = np.array([DT], dtype=float)
        else:
            DT = np.asarray(DT, dtype=float)

        if isinstance(Omega, (float, int)):
            Omega = np.array([Omega], dtype=float)
        else:
            Omega = np.asarray(Omega, dtype=float)

        self.DT: npt.NDArray[np.float64] = DT
        "Temperature difference between the plates [K]"

        self.Omega: npt.NDArray[np.float64] = Omega
        "Rotation rate along the cell axis [rad/s]"

        self.Pr = self.fluid.Pr
        """Prandtl number"""

        self.Ra = (
            DT * fluid.alpha * cell.g * cell.L**3 / (fluid.kappa * fluid.nu)
        )
        with np.errstate(divide="ignore"):
            # Numpy handles divide by 0 alright by setting it to `np.inf`, so
            # it's okay to ignore the "RuntimeWarning: divide by zero" and
            # silence it.
            """Rayleigh number"""
            self.Ek = fluid.nu / (2 * Omega * cell.L**2)
            """Ekman number"""
            self.Ro = np.sqrt(cell.g * fluid.alpha * DT / cell.L) / (2 * Omega)
            """Rossby number"""
            self.Ta = 4 * cell.L**4 * Omega**2 / fluid.nu**2
            """Taylor number"""
            self.Fr = Omega**2 * (cell.D / 2) / cell.g
            """Centrifugal Froude number"""
            self.Ra_c = 7.8 * np.power(self.Ek, -4 / 3)
            """Onset of convection with fixed boundaries: Ra_c = 7.8 * Ek^(-4/3)
            """

    def report(self):
        """Print all fluid properties, the RBC cell, the driving parameters and
        their non-dimensional numbers as a table to the terminal. Will only
        print the first element if arrays were used as input.
        """
        self.fluid.report()
        print(f"  RBC cell: {self.cell.name}")
        print(
            f"  └ D = {self.cell.D:5.3f} m, L = {self.cell.L:<5.3f} m, "
            f"g = {self.cell.g:<4.2f} m/s^2"
        )
        print(f"  @ Temperature diff. | DT    = {self.DT[0]:<6.3f} K")
        print(
            f"  @ Rotation rate     | Omega = {self.Omega[0]:<6.4f} rad/s = "
            f"{self.Omega[0] / np.pi / 2 * 60:6.3f} rpm"
        )
        print(fluidprop.HRULE)
        pprint("Rayleigh", "Ra", self.Ra[0])

        if not np.all(self.Omega == 0):
            # fmt: off
            pprint("Ekman"      , "Ek"  , self.Ek[0]  , format_spec="e")
            pprint("Rossby"     , "Ro"  , self.Ro[0]  , format_spec="f")
            with np.errstate(divide="ignore"):
                pprint("Inv. Rossby", "1/Ro", 1/self.Ro[0], format_spec="f")
            pprint("Taylor"     , "Ta"  , self.Ta[0]  , format_spec="e")
            pprint("Centrifugal Froude", "Fr"  , self.Fr[0]  , format_spec="f")
            print("    Onset convect. fixed   | Ra_c    = 7.8 * Ek^(-4/3)")
            pprint(
                "Onset convect. fixed", "Ra_c"   , self.Ra_c[0],
                format_spec="e"
            )
            pprint(
                "N times from onset"  , "Ra/Ra_c", self.Ra[0] / self.Ra_c[0],
                format_spec="f", N_decimals=1
            )
            # fmt: on

        print(fluidprop.HRULE)


# ------------------------------------------------------------------------------
#   show_cli
# ------------------------------------------------------------------------------


def show_cli() -> RBC:
    """Show a minimal command-line interface which guides the user to enter all
    relevant parameters for Rayleigh-Bénard Convection.

    Returns:
        The resulting `RBC()` class instance.
    """
    # Table of all known RBC cells
    RBC_cell_table = (
        "   # | name                 |   D   |   L   | Gamma| g\n"
        "     |                      |   m   |   m   |   -  | m/s^2\n"
        "  ----------------------------------------------------------\n"
    )
    for cell_idx, cell in enumerate(RBC_cells):
        RBC_cell_table += (
            f"  {cell_idx:>2} | {cell.name:<20} | {cell.D:<5.3f} | "
            f"{cell.L:<5.3f} | {cell.Gamma:<4.2f} | {cell.g:<4.2f}\n"
        )

    DT = np.nan
    """Temperature difference between the plates [K]"""
    Omega = np.nan
    """Rotation rate along the cell axis [rad/s]"""

    # ---------------------------------
    #   Ask user for input parameters
    # ---------------------------------

    fluid = fluidprop.show_cli()

    while True:
        # -------------------
        #   Select RBC cell
        # -------------------

        print("\nKnown RBC cells:")
        print(RBC_cell_table)

        ans = input("Enter cell number: ")
        try:
            RBC_cell_idx = int(ans)
        except ValueError:
            print("Not a valid number.")
            continue

        valid_values = range(len(RBC_cells))
        if RBC_cell_idx not in valid_values:
            print("Not a valid number.")
            continue

        cell = RBC_cells[RBC_cell_idx]

        # --------------------------------
        #   Enter temperature difference
        # --------------------------------

        ans = input("Enter temperature diff. | DT    [K]    : ")
        if ans.lower() == "f":
            ans = input("Enter temperature diff. | DT    ['F]   : ")
            try:
                DT = (float(ans) - 32) * 5 / 9
            except ValueError:
                print("Not a valid number.")
                continue

        else:
            try:
                DT = float(ans)
            except ValueError:
                print("Not a valid number.")
                continue

        # -----------------------
        #   Enter rotation rate
        # -----------------------

        ans = input("Enter rotation rate     | Omega [rad/s]: ")
        if ans.lower() == "r":
            ans = input("Enter rotation rate     | Omega [rpm]  : ")
            try:
                Omega = float(ans) / 60 * 2 * np.pi
            except ValueError:
                print("Not a valid number. Rotation rate set to 0.")
                Omega = 0
                break

        elif ans.lower() in ("h", "v"):  # [Hz] or [rev/s]
            ans = input("Enter rotation rate     | Omega [rev/s]: ")
            try:
                Omega = float(ans) * 2 * np.pi
            except ValueError:
                print("Not a valid number. Rotation rate set to 0.")
                Omega = 0
                break

        else:
            try:
                Omega = float(ans)
            except ValueError:
                print("Not a valid number. Rotation rate set to 0.")
                Omega = 0
                break

        # Made it successfully till the end
        break

    return RBC(fluid=fluid, cell=cell, DT=DT, Omega=Omega)


# ------------------------------------------------------------------------------
#   main
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    rbc = show_cli()
    print()
    rbc.report()
