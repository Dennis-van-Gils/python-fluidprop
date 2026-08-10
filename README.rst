|pypi| |python| |black| |license|

.. |pypi| image:: https://img.shields.io/pypi/v/fluidprop
    :target: https://pypi.org/project/fluidprop
.. |python| image:: https://img.shields.io/pypi/pyversions/fluidprop
    :target: https://pypi.org/project/fluidprop
.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
.. |license| image:: https://img.shields.io/badge/License-MIT-purple.svg
    :target: https://github.com/Dennis-van-Gils/python-dvg-devices/blob/master/LICENSE.txt

fluidprop
=========
Easy access to thermodynamic fluid properties as a function of temperature and
pressure. Comes with a minimal command-line interface for quick inspection.
Provides class ``FluidProperties()`` useful for working out dataseries in your
own scripts.

- Github: https://github.com/Dennis-van-Gils/python-fluidprop
- PyPI: https://pypi.org/project/fluidprop

Installation::

    pip install fluidprop

or if you're on macOS or Linux, try::

    pip3 install fluidprop

Raison d'être
-------------

Science advances every year and with it the accuracy of tabulated/parametrized
fluid properties. Don't reinvent the wheel by coding your own equation-of-state
models from the literature. Or copy-pasting possibly outdated values from the
internet.

EOS models
----------

Thermodynamic properties are provided by CoolProp, the open-source alternative
to `NIST refprop <https://www.nist.gov/srd/refprop>`_, with most of the
calculations relying on the same equation-of-state (EOS) models as refprop.

* http://www.coolprop.org/
* http://pubs.acs.org/doi/abs/10.1021/ie4033999

Command-line interface
======================

You can run this module from the terminal with::

    python -m fluidprop

or if you're on macOS or Linux, try::

    python3 -m fluidprop

It will show a minimal command-line interface which guides the user to enter a
fluid, temperature and pressure. It will print out its thermodynamic properties
as a table to the terminal.

BONUS for the Rayleigh-Bénard convection community. Running::

    python -m rbc

will show a command-line interface acting as a 'pocket calculator' to calculate
the Rayleigh and other numbers based on the user input.

Example output of `fluidprop`::

    https://github.com/Dennis-van-Gils/python-fluidprop
    Thermodynamic properties by CoolProp v6.6.0
    http://pubs.acs.org/doi/abs/10.1021/ie4033999

    All known fluids:
        0 | 1-Butene             46 | MethylLinoleate      92 | R1234yf
        1 | Acetone              47 | MethylLinolenate     93 | R1234ze(E)
        2 | Air                  48 | MethylOleate         94 | R1234ze(Z)
        3 | Ammonia              49 | MethylPalmitate      95 | R124
        4 | Argon                50 | MethylStearate       96 | R1243zf
        5 | Benzene              51 | MM                   97 | R125
        6 | CarbonDioxide        52 | n-Butane             98 | R13
        7 | CarbonMonoxide       53 | n-Decane             99 | R1336mzz(E)
        8 | CarbonylSulfide      54 | n-Dodecane          100 | R1336mzz(Z)
        9 | Chlorine             55 | n-Heptane           101 | R134a
       10 | cis-2-Butene         56 | n-Hexane            102 | R13I1
       11 | CycloHexane          57 | n-Nonane            103 | R14
       12 | Cyclopentane         58 | n-Octane            104 | R141b
       13 | CycloPropane         59 | n-Pentane           105 | R142b
       14 | D4                   60 | n-Perfluorobutane   106 | R143a
       15 | D5                   61 | n-Perfluorohexane   107 | R152A
       16 | D6                   62 | n-Perfluoropentane  108 | R161
       17 | Deuterium            63 | n-Propane           109 | R21
       18 | Dichloroethane       64 | n-Undecane          110 | R218
       19 | DiethylEther         65 | Neon                111 | R22
       20 | DimethylCarbonate    66 | Neopentane          112 | R227EA
       21 | DimethylEther        67 | Nitrogen            113 | R23
       22 | Ethane               68 | NitrousOxide        114 | R236EA
       23 | Ethanol              69 | Novec649            115 | R236FA
       24 | EthylBenzene         70 | o-Xylene            116 | R245ca
       25 | Ethylene             71 | OrthoDeuterium      117 | R245fa
       26 | EthyleneOxide        72 | OrthoHydrogen       118 | R32
       27 | Fluorine             73 | Oxygen              119 | R365MFC
       28 | HeavyWater           74 | p-Xylene            120 | R40
       29 | Helium               75 | ParaDeuterium       121 | R404A
       30 | HFE143m              76 | ParaHydrogen        122 | R407C
       31 | Hydrogen             77 | Propylene           123 | R41
       32 | HydrogenChloride     78 | PropyleneGlycol     124 | R410A
       33 | HydrogenSulfide      79 | Propyne             125 | R507A
       34 | IsoButane            80 | R11                 126 | RC318
       35 | IsoButene            81 | R1123               127 | SES36
       36 | Isohexane            82 | R113                128 | SulfurDioxide
       37 | Isopentane           83 | R1130(E)            129 | SulfurHexafluoride
       38 | Krypton              84 | R1132(E)            130 | Tetrahydrofuran
       39 | m-Xylene             85 | R114                131 | Toluene
       40 | MD2M                 86 | R115                132 | trans-2-Butene
       41 | MD3M                 87 | R116                133 | VinylChloride
       42 | MD4M                 88 | R12                 134 | Water
       43 | MDM                  89 | R1224YDZ            135 | Xenon
       44 | Methane              90 | R123
       45 | Methanol             91 | R1233zd(E)

    Enter fluid number: 121
    Enter temperature | T ['C]  : 20
    Enter pressure    | P [bar] : a
    Enter pressure    | P [atm] : 1

    ------------------------------------------------------------
    Liquid: Water (H₂O)
    @ Temperature | T =   20.000 'C = 293.150 K
    @ Pressure    | P =    1.013 bar
    ------------------------------------------------------------
        Molecular weight       | MW      = 18.01527    g/mol
        Density                | rho     = 9.982e+02   kg/m^3
        Kinematic viscosity    | nu      = 1.003e-06   m^2/s
        Dynamic   viscosity    | eta     = 1.002e-03   kg/(m s)
        Thermal exp. coeff.    | alpha   = 2.068e-04   1/K
        Thermal diffusivity    | kappa   = 1.432e-07   m^2/s
        Thermal conductivity   | lambda_ = 5.980e-01   W/(m K)
        Isobaric  heat capac.  | Cp      = 4.184e+03   J/(kg K)
        Isochoric heat capac.  | Cv      = 4.157e+03   J/(kg K)
        Isothermal compress.   | comp    = 4.589e-10   1/Pa
        Prandtl                | Pr      = 7.008
    ------------------------------------------------------------

When asked to enter the temperature in ``['C]``, you can *once* enter a single
character instead to change the input unit to::

    k | [K]     Kelvin                  K - 273.15 'C
    f | ['F]    Degrees Fahrenheit      ('F - 32) * 5 / 9 'C

When asked to enter the pressure in ``[bar]``, you can *once* enter a single
character instead to change the input unit to::

    a | [atm]   Atmosphere              = 1.01325 bar
    m | [mmHg]  Millimeter mercury      ≈ 1 atm / 760
    p | [psi]   Pounds per square inch  = 1 / 14.504 bar
    t | [torr]  Torr                    = 1 atm / 760

FluidProperties()
=================

This class evaluates thermodynamic fluid properties of the given fluid at the
given temperature(s) in ``['C]`` and pressure(s) in ``[bar]``. The results are
stored as properties to this class as ``numpy.ndarray`` arrays. Useful for
working out dataseries.

Example:

.. code-block:: python

    from fluidprop import FluidProperties

    fluid = FluidProperties("Water", 20, 1)
    print(fluid.rho)  # [998.2065435]

    fluid = FluidProperties("Water", [20, 21, 22], 1)
    print(fluid.rho)  # [998.2065435  997.99487638 997.77288644]

List of stored properties::

    coolprop_name (str): CoolProp name of the fluid.

    formula       (str): Chemical formula of the fluid.

    MW      (float)  : Molecular weight               [kg/mol]

    T       (ndarray): Evaluated temperature          [K]

    P       (ndarray): Evaluated pressure             [Pa]

    rho     (ndarray): Density                        [kg/m^3]

    nu      (ndarray): Kinematic viscosity            [m^2/s]

    eta     (ndarray): Dynamic/shear viscosity        [kg/(m s)]

    alpha   (ndarray): Thermal expansion coefficient  [1/K]

    kappa   (ndarray): Thermal diffusivity            [m^2/s]

    lambda_ (ndarray): Thermal conductivity           [W/(m K)]

    Cp      (ndarray): Isobaric heat capacity         [J/(kg K)]

    Cv      (ndarray): Isochoric heat capacity        [J/(kg K)]

    comp    (ndarray): Isothermal compressibility     [1/Pa]

    Pr      (ndarray): Prandtl number                 [-]

Dennis van Gils, 10-08-2026
