#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Easy access to thermodynamic fluid properties as a function of temperature
and pressure. Comes with a minimal command-line interface for quick inspection.
Provides class ``FluidProperties()`` useful for working out dataseries in your
own scripts.

Thermodynamic properties are provided by CoolProp:
* http://www.coolprop.org/
* http://pubs.acs.org/doi/abs/10.1021/ie4033999
"""
__author__ = "Dennis van Gils"
__authoremail__ = "vangils.dennis@gmail.com"
__url__ = "https://github.com/Dennis-van-Gils/python-fluidprop"
__date__ = "12-05-2024"
__version__ = "1.1.0"

import re
from typing import Union

import numpy.typing as npt
import numpy as np
import CoolProp
import CoolProp.CoolProp as CP


ZERO_C = 273.15
"""0 'C in K"""
P_ATM = 1.01325
"""bar per 1 atm"""
P_PSI = 1 / 14.504
"""bar per 1 psi"""
HRULE = "-" * 60
"""Horizontal rule"""

# fmt: off
FLUID_SELECTION = [
    ("Air"                  , "mixture"),
    ("Hydrogen"             , "H_{2}"),
    ("Helium"               , "He"),
    ("Nitrogen"             , "N_{2}"),
    ("Oxygen"               , "O_{2}"),
    ("CarbonDioxide"        , "CO_{2}"),
    ("SulfurHexafluoride"   , "SF_{6}"),
    ("Water"                , "H_{2}O"),
    ("HeavyWater"           , "D_{2}O"),
    ("Methanol"             , "CH_{3}OH"),
    ("Ethanol"              , "C_{2}H_{5}OH"),
    # ("Acetone"              , "C_{3}H_{6}O"),
]
# fmt: on
"""Subselection of the most common fluids to chose from.
It is a list of tuples, where each tuple contains 2 strings as follows ::

    tuple[0] (str): CoolProp name
    tuple[1] (str): Chemical name
"""

# ------------------------------------------------------------------------------
#   FluidProperties
# ------------------------------------------------------------------------------


class FluidProperties:
    """Evaluates thermodynamic fluid properties of the given fluid at the given
    temperature(s) in ``['C]`` and pressure(s) in ``[bar]``. The results are
    stored as properties to this class as ``numpy.ndarray`` arrays. Useful for
    working out dataseries.

    Example ::

        fluid = FluidProperties("Water", 20, 1)
        print(fluid.rho)  # [998.2065435]

        fluid = FluidProperties("Water", [20, 21, 22], 1)
        print(fluid.rho)  # [998.2065435  997.99487638 997.77288644]

    Args:
        coolprop_name (`str`):
            The CoolProp name of the fluid to evaluate.

        T_in_deg_C (`float` | `list[float]` | `numpy.ndarray[]`):
            Temperature ['C] to evaluate fluid properties at.

        P_in_bar (`float` | `list[float]` | `numpy.ndarray[]`):
            Pressure [bar] to evaluate fluid properties at.

    Properties:
        coolprop_name (`str`):
            CoolProp name of the fluid.

        formula (`str`):
            Chemical formula of the fluid.

        MW (`float`):
            Molecular weight [kg/mol].

        T (`numpy.ndarray[]`):
            Evaluated temperature [K].

        P (`numpy.ndarray[]`):
            Evaluated pressure [Pa].

        rho (`numpy.ndarray[]`):
            Density [kg/m^3].

        nu (`numpy.ndarray[]`):
            Kinematic viscosity [m^2/s].

        eta (`numpy.ndarray[]`):
            Dynamic/shear viscosity [kg/(m s)].

        alpha (`numpy.ndarray[]`):
            Thermal expansion coefficient [1/K].

        kappa (`numpy.ndarray[]`):
            Thermal diffusivity [m^2/s].

        lambda_ (`numpy.ndarray[]`):
            Thermal conductivity [W/(m K)].

        Cp (`numpy.ndarray[]`):
            Isobaric heat capacity [J/(kg K)].

        Cv (`numpy.ndarray[]`):
            Isochoric heat capacity [J/(kg K)].

        comp (`numpy.ndarray[]`):
            Isothermal compressibility [1/Pa].

        Pr (`numpy.ndarray[]`):
            Prandtl number.
    """

    def __init__(
        self,
        coolprop_name: str,
        T_in_deg_C: Union[float, list[float], npt.NDArray],
        P_in_bar: Union[float, list[float], npt.NDArray],
    ):
        # -------------------------
        #   Check input arguments
        # -------------------------

        if isinstance(T_in_deg_C, (float, int)):
            T_in_deg_C = np.array([T_in_deg_C], dtype=float)
        else:
            T_in_deg_C = np.asarray(T_in_deg_C, dtype=float)

        if isinstance(P_in_bar, (float, int)):
            P_in_bar = np.array([P_in_bar], dtype=float)
        else:
            P_in_bar = np.asarray(P_in_bar, dtype=float)

        if T_in_deg_C.ndim != 1 or P_in_bar.ndim != 1:
            raise ValueError(
                "Arguments `T_in_deg_C` and `P_in_bar` must be one-dimensional "
                "arrays."
            )

        if len(T_in_deg_C) > 1 and len(P_in_bar) == 1:
            P_in_bar = np.repeat(P_in_bar, len(T_in_deg_C))

        if len(P_in_bar) > 1 and len(T_in_deg_C) == 1:
            T_in_deg_C = np.repeat(T_in_deg_C, len(P_in_bar))

        if len(T_in_deg_C) != len(P_in_bar):
            raise ValueError(
                "Arguments `T_in_deg_C` and `P_in_bar` have unequal lengths."
            )

        array_len = np.size(T_in_deg_C)
        nan_array = np.empty(array_len)
        nan_array[:] = np.nan

        # -----------
        #   Members
        # -----------

        self.coolprop_name: str = coolprop_name
        """CoolProp name of the fluid"""

        self.formula: str = ""
        """Chemical formula of the fluid"""

        self.MW: float = np.nan
        """Molecular weight [kg/mol]"""

        self.T: npt.NDArray[np.float64] = np.add(T_in_deg_C, ZERO_C)
        """Evaluated temperature [K]"""

        self.P: npt.NDArray[np.float64] = np.multiply(P_in_bar, 1e5)
        """Evaluated pressure [Pa]"""

        self.rho: npt.NDArray[np.float64] = np.copy(nan_array)
        """Density [kg/m^3]"""

        self.nu: npt.NDArray[np.float64] = np.copy(nan_array)
        """Kinematic viscosity [m^2/s]"""

        self.eta: npt.NDArray[np.float64] = np.copy(nan_array)
        """Dynamic/shear viscosity [kg/(m s)]"""

        self.alpha: npt.NDArray[np.float64] = np.copy(nan_array)
        """Thermal expansion coefficient [1/K]"""

        self.kappa: npt.NDArray[np.float64] = np.copy(nan_array)
        """Thermal diffusivity [m^2/s]"""

        self.lambda_: npt.NDArray[np.float64] = np.copy(nan_array)
        """Thermal conductivity [W/(m K)]"""

        self.Cp: npt.NDArray[np.float64] = np.copy(nan_array)
        """Isobaric heat capacity [J/(kg K)]"""

        self.Cv: npt.NDArray[np.float64] = np.copy(nan_array)
        """Isochoric heat capacity [J/(kg K)]"""

        self.comp: npt.NDArray[np.float64] = np.copy(nan_array)
        """Isothermal compressibility [1/Pa]"""

        self.Pr: npt.NDArray[np.float64] = np.copy(nan_array)
        """Prandtl number"""

        # ------------------------------
        #   Calculate fluid properties
        # ------------------------------

        found_match = False
        for item_ in FLUID_SELECTION:
            if item_[0] == coolprop_name:
                found_match = True
                self.formula = item_[1]
                break

        if not found_match:
            self.formula = CP.get_fluid_param_string(coolprop_name, "formula")
            self.formula = self.formula.replace("_{1}", "")

        # Molecular weight [kg/mol]
        self.MW = CP.PropsSI(coolprop_name, "M")

        # http://coolprop.org/coolprop/HighLevelAPI.html#table-of-string-inputs-to-propssi-function
        requested_quantities = [
            "DMASS",
            "VISCOSITY",
            "ISOBARIC_EXPANSION_COEFFICIENT",
            "CONDUCTIVITY",
            "CPMASS",
            "CVMASS",
            "ISOTHERMAL_COMPRESSIBILITY",
            "PRANDTL",
        ]

        for idx_, (T, P) in enumerate(zip(self.T, self.P)):
            values = np.zeros(len(requested_quantities))
            values[:] = np.nan

            for quantity_idx, quantity in enumerate(requested_quantities):
                try:
                    value = CP.PropsSI(quantity, "T", T, "P", P, coolprop_name)
                except ValueError as e:
                    value = np.nan
                    print(e)

                values[quantity_idx] = value

            # fmt: off
            (
                self.rho[idx_],     # Density                       [kg/m^3]
                self.eta[idx_],     # Dynamic viscosity             [kg/(m s)]
                self.alpha[idx_],   # Thermal expansion coefficient [1/K]
                self.lambda_[idx_], # Thermal conductivity          [W/(m K)]
                self.Cp[idx_],      # Isobaric heat capacity        [J/(kg K)]
                self.Cv[idx_],      # Isochoric heat capacity       [J/(kg K)]
                self.comp[idx_],    # Isothermal compressibility    [1/Pa]
                self.Pr[idx_],      # Prandtl number                [-]
            ) = values
            # fmt: on

        # Derived: Kinematic viscosity [m^2/s]
        self.nu = self.eta / self.rho
        # Derived: Thermal diffusivity [m^2/s]
        self.kappa = self.lambda_ / self.rho / self.Cp

    def report(self):
        """Print all fluid properties as a table to the terminal. Will only
        print the first element if arrays were used as input.
        """

        print(HRULE)
        print(
            f"  Liquid: {self.coolprop_name} ({utf8_subscripts(self.formula)})"
        )
        print(
            f"  @ Temperature | T = {self.T[0] - ZERO_C:8.3f} 'C = "
            f"{self.T[0]:.3f} K"
        )
        print(f"  @ Pressure    | P = {self.P[0]/1e5:8.3f} bar")
        print(HRULE)

        # fmt: off
        pprint("Molecular weight"      , "MW"     , self.MW * 1e3  , "g/mol",
               format_spec="f", N_decimals=5)
        pprint("Density"               , "rho"    , self.rho[0]    , "kg/m^3")
        pprint("Kinematic viscosity"   , "nu"     , self.nu[0]     , "m^2/s")
        pprint("Dynamic   viscosity"   , "eta"    , self.eta[0]    , "kg/(m s)")
        pprint("Thermal exp. coeff."   , "alpha"  , self.alpha[0]  , "1/K")
        pprint("Thermal diffusivity"   , "kappa"  , self.kappa[0]  , "m^2/s")
        pprint("Thermal conductivity"  , "lambda_", self.lambda_[0], "W/(m K)")
        pprint("Isobaric  heat capac." , "Cp"     , self.Cp[0]     , "J/(kg K)")
        pprint("Isochoric heat capac." , "Cv"     , self.Cv[0]     , "J/(kg K)")
        pprint("Isothermal compress. " , "comp"   , self.comp[0]   , "1/Pa")
        pprint("Prandtl"               , "Pr"     , self.Pr[0]     ,
               format_spec="f")
        # fmt: on
        print(HRULE)


# ------------------------------------------------------------------------------
#   Helper functions
# ------------------------------------------------------------------------------


def pprint(
    descr: str,
    abbrev: str,
    value: float,
    unit: str = "",
    format_spec: str = "e",
    N_decimals: int = 3,
):
    """Pretty print in columns."""
    print(
        f"    {descr:<22} | {abbrev:<7} = "
        f"{value:<11.{N_decimals}{format_spec}} {unit:<8}"
    )


def utf8_subscripts(text: str) -> str:
    """Find all instances of '_{###}' and replace them with UTF-8 subscript
    numbers, except for '_{1}' which we will remove."""
    subscript_dict = {
        "0": "\u2080",  # ₀
        "1": "\u2081",  # ₁
        "2": "\u2082",  # ₂
        "3": "\u2083",  # ₃
        "4": "\u2084",  # ₄
        "5": "\u2085",  # ₅
        "6": "\u2086",  # ₆
        "7": "\u2087",  # ₇
        "8": "\u2088",  # ₈
        "9": "\u2089",  # ₉
    }

    text = text.replace("_{1}", "")
    text = re.sub(
        r"_{(\d+)}",
        lambda match: "".join(subscript_dict[i] for i in match.group(1)),
        text,
    )

    return text


# ------------------------------------------------------------------------------
#   main
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    print(__url__)
    print(f"Thermodynamic properties by CoolProp v{CoolProp.__version__}")
    print("http://pubs.acs.org/doi/abs/10.1021/ie4033999")

    # Table of predefined subselection of fluids
    table_1 = ""
    for idx, item in enumerate(FLUID_SELECTION):
        table_1 += f"  {idx:>3} | {utf8_subscripts(item[1]):<8} | {item[0]}\n"

    # Table of all known fluids
    table_2 = ""
    N_fluids = len(CoolProp.__fluids__)
    N_rows = int(np.ceil(N_fluids / 3))  # Spread out over 3 columns
    for idx_1 in range(N_rows):
        idx_2 = idx_1 + N_rows
        idx_3 = idx_2 + N_rows

        table_2 += f"  {idx_1:>3} | {CoolProp.__fluids__[idx_1]:<18}"
        if idx_2 < N_fluids:
            table_2 += f"  {idx_2:>3} | {CoolProp.__fluids__[idx_2]:<18}"
        if idx_3 < N_fluids:
            table_2 += f"  {idx_3:>3} | {CoolProp.__fluids__[idx_3]:<18}"
        table_2 += "\n"

    # ---------------------------------
    #   Ask user for input parameters
    # ---------------------------------

    while True:
        # ----------------
        #   Select fluid
        # ----------------

        print("\nQuick selection fluids:")
        print(table_1)
        print("Enter nothing to show all fluids.")

        ans = input("Enter fluid number: ")
        if ans == "":
            print("All known fluids:")
            print(table_2)

            ans = input("Enter fluid number: ")
            try:
                fluid_idx = int(ans)
            except ValueError:
                fluid_idx = -1

            if fluid_idx not in range(len(CoolProp.__fluids__)):
                print("Not a valid number.")
                continue

            fluid_name = CoolProp.__fluids__[fluid_idx]

        else:
            try:
                fluid_idx = int(ans)
            except ValueError:
                fluid_idx = -1

            if fluid_idx not in range(len(FLUID_SELECTION)):
                print("Not a valid number.")
                continue

            fluid_name = FLUID_SELECTION[fluid_idx][0]

        # ---------------------
        #   Enter temperature
        # ---------------------

        ans = input("Enter temperature | T ['C]  : ")
        if ans.lower() == "f":
            ans = input("Enter temperature | T ['F]  : ")
            try:
                temperature_deg_C = (float(ans) - 32) * 5 / 9
            except ValueError:
                print("Not a valid number.")
                continue

        elif ans.lower() == "k":
            ans = input("Enter temperature | T [K]   : ")
            try:
                temperature_deg_C = float(ans) - ZERO_C
            except ValueError:
                print("Not a valid number.")
                continue

        else:
            try:
                temperature_deg_C = float(ans)
            except ValueError:
                print("Not a valid number.")
                continue

        # ------------------
        #   Enter pressure
        # ------------------

        ans = input("Enter pressure    | P [bar] : ")
        if ans.lower() == "a":
            ans = input("Enter pressure    | P [atm] : ")
            try:
                pressure_bar = float(ans) * P_ATM
            except ValueError:
                print("Not a valid number.")
                continue

        elif ans.lower() == "m":
            ans = input("Enter pressure    | P [mmHg]: ")
            try:
                pressure_bar = float(ans) * P_ATM / 760
            except ValueError:
                print("Not a valid number.")
                continue

        elif ans.lower() == "p":
            ans = input("Enter pressure    | P [psi] : ")
            try:
                pressure_bar = float(ans) * P_PSI
            except ValueError:
                print("Not a valid number.")
                continue

        elif ans.lower() == "t":
            ans = input("Enter pressure    | P [torr]: ")
            try:
                pressure_bar = float(ans) * P_ATM / 760
            except ValueError:
                print("Not a valid number.")
                continue

        else:
            try:
                pressure_bar = float(ans)
            except ValueError:
                print("Not a valid number.")
                continue

        # Made it successfully till the end
        break

    # ----------------------------------
    #   Calculate and report to screen
    # ----------------------------------

    fluid = FluidProperties(
        coolprop_name=fluid_name,
        T_in_deg_C=temperature_deg_C,
        P_in_bar=pressure_bar,
    )

    print()
    fluid.report()
