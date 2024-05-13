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
        0 | 1-Butene             41 | MD4M                 82 | R1233zd(E)
        1 | Acetone              42 | MDM                  83 | R1234yf
        2 | Air                  43 | Methane              84 | R1234ze(E)
        3 | Ammonia              44 | Methanol             85 | R1234ze(Z)
        4 | Argon                45 | MethylLinoleate      86 | R124
        5 | Benzene              46 | MethylLinolenate     87 | R1243zf
        6 | CarbonDioxide        47 | MethylOleate         88 | R125
        7 | CarbonMonoxide       48 | MethylPalmitate      89 | R13
        8 | CarbonylSulfide      49 | MethylStearate       90 | R134a
        9 | cis-2-Butene         50 | MM                   91 | R13I1
       10 | CycloHexane          51 | n-Butane             92 | R14
       11 | Cyclopentane         52 | n-Decane             93 | R141b
       12 | CycloPropane         53 | n-Dodecane           94 | R142b
       13 | D4                   54 | n-Heptane            95 | R143a
       14 | D5                   55 | n-Hexane             96 | R152A
       15 | D6                   56 | n-Nonane             97 | R161
       16 | Deuterium            57 | n-Octane             98 | R21
       17 | Dichloroethane       58 | n-Pentane            99 | R218
       18 | DiethylEther         59 | n-Propane           100 | R22
       19 | DimethylCarbonate    60 | n-Undecane          101 | R227EA
       20 | DimethylEther        61 | Neon                102 | R23
       21 | Ethane               62 | Neopentane          103 | R236EA
       22 | Ethanol              63 | Nitrogen            104 | R236FA
       23 | EthylBenzene         64 | NitrousOxide        105 | R245ca
       24 | Ethylene             65 | Novec649            106 | R245fa
       25 | EthyleneOxide        66 | o-Xylene            107 | R32
       26 | Fluorine             67 | OrthoDeuterium      108 | R365MFC
       27 | HeavyWater           68 | OrthoHydrogen       109 | R40
       28 | Helium               69 | Oxygen              110 | R404A
       29 | HFE143m              70 | p-Xylene            111 | R407C
       30 | Hydrogen             71 | ParaDeuterium       112 | R41
       31 | HydrogenChloride     72 | ParaHydrogen        113 | R410A
       32 | HydrogenSulfide      73 | Propylene           114 | R507A
       33 | IsoButane            74 | Propyne             115 | RC318
       34 | IsoButene            75 | R11                 116 | SES36
       35 | Isohexane            76 | R113                117 | SulfurDioxide
       36 | Isopentane           77 | R114                118 | SulfurHexafluoride
       37 | Krypton              78 | R115                119 | Toluene
       38 | m-Xylene             79 | R116                120 | trans-2-Butene
       39 | MD2M                 80 | R12                 121 | Water
       40 | MD3M                 81 | R123                122 | Xenon

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

Dennis van Gils, 13-05-2024
