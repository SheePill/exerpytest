import logging

import numpy as np

from exerpy.components.component import Component
from exerpy.components.component import component_registry


@component_registry
class SolarTower(Component):
    r"""
    class for exergy and exergoeconomic analysis of Solar Tower.

    This class performs exergy and exergoeconomic analysis calculations for Solar Tower components,
    accounting for one LOGIC, one inlet and one outlet stream across various temperature regimes, including
    above ambient temperature, and optional dissipative behavior.

    Attributes
    ----------
    E_F : float
        Exergy fuel of the component :math:`\dot{E}_\mathrm{F}` in :math:`\mathrm{W}`.
    E_P : float
        Exergy product of the component :math:`\dot{E}_\mathrm{P}` in :math:`\mathrm{W}`.
    E_D : float
        Exergy destruction of the component :math:`\dot{E}_\mathrm{D}` in :math:`\mathrm{W}`.
    epsilon : float
        Exergetic efficiency of the component :math:`\varepsilon` in :math:`-`.
    inl : dict
        Dictionary containing inlet stream data with mass flows and specific exergies.
    outl : dict
        Dictionary containing outlet stream data with mass flows and specific exergies.
    logic_inl  : dict
        Dictionary containing logic connection to the heliostat field.
    """

    def __init__(self, **kwargs):
        r"""
        Initializes the SolarTower component.
        """
        super().__init__(**kwargs)

    def calc_exergy_balance(self, T0: float, p0: float, split_physical_exergy) -> None:
        r"""
        Method to calculate the exergy balance of the Solar Tower component.

        Parameters
        ----------
        T0 : float
            Ambient temperature in :math:`\mathrm{K}`.
        p0 : float
            Ambient pressure in :math:`\mathrm{Pa}`.
        split_physical_exergy : bool
            Flag to indicate whether to split physical exergy into thermal and mechanical parts.

        Returns
        -------
        None
            The method updates the attributes E_F, E_P, E_D, and epsilon of the component instance.

        Raises
        ------
        ValueError
            If the number of inlet or outlet streams is not equal to one.

        Notes
        -----
        The exergy balance is calculated using the following equations:

        - Exergy fuel: :math:`\dot{E}_\mathrm{F} = \sum \dot{m}_\mathrm{in} \cdot e_\mathrm{in}`
        - Exergy product: :math:`\dot{E}_\mathrm{P} = \sum \dot{m}_\mathrm{out} \cdot e_\mathrm{out}`
        - Exergy destruction: :math:`\dot{E}_\mathrm{D} = \dot{E}_\mathrm{F} - \dot{E}_\mathrm{P}`
        - Exergetic efficiency: :math:`\varepsilon = \frac{\dot{E}_\mathrm{P}}{\dot{E}_\mathrm{F}}`

        where :math:`\dot{m}` is the mass flow rate and :math:`e` is the specific exergy of the streams.
        """
        # Validate the number of inlets and outlets
        if not hasattr(self, 'inl') or not hasattr(self, 'outl') or not hasattr(self, 'logic_inl') or len(self.inl) != 1 or len(self.outl) != 1 or len(self.logic_inl) != 1:
            msg = "SolarTower requires exactly one inlet, one outlet, and one logic inlet to Heliostatfield."
            logging.error(msg)
            raise ValueError(msg)
    
        # Extract inlet, outlet streams and one logic inlet
        inlet = self.inl[0]
        outlet = self.outl[0]
        logic_inlet = self.logic_inl[0]
    
        # For solar tower, Q comes from the logic inlet (solar energy input)
        #Q = logic_inlet['Q']  # Solar heat input from heliostat field
        #Q = outlet['m'] * outlet['h'] - inlet['m'] * inlet['h']

        # Initialize E_P and E_F
        self.E_P = 0.0
        self.E_F = 0.0

        # Case 1: Heat is added to the system (Q > 0)
        if Q > 0:
            if inlet['T'] >= T0 and outlet['T'] >= T0:
                if split_physical_exergy:
                    self.E_P = outlet['m'] * (outlet['e_PH'] - inlet['e_PH'])
                    self.E_F = logic_inlet['Q']  # Solar heat input
                else:
                    self.E_P = outlet['m'] * (outlet['e_PH'] - inlet['e_PH'])
                    self.E_F = logic_inlet['Q']  
            
            elif inlet['T'] < T0 and outlet['T'] > T0:
                if split_physical_exergy:
                    self.E_P = outlet['m'] * (outlet['e_T'] + inlet['e_T'])
                    self.E_F = logic_inlet['Q']
                else:
                    self.E_P = outlet['m'] * (outlet['e_PH'] - inlet['e_PH'])
                    self.E_F = logic_inlet['Q']

            elif inlet['T'] < T0 and outlet['T'] < T0:
                logging.warning(
                    "The inlet and outlet temperatures are both below ambient (T_in < T0 and T_out < T0). This case is not typically expected for a Solar Tower."
                )
                self.E_P = np.nan
                self.E_F = np.nan
            else:
                logging.warning(
                    "SolarTower: unimplemented case (Q > 0, T_in > T0 > T_out?)."
                )
                self.E_P = np.nan
                self.E_F = outlet['m'] * (inlet['e_PH'] - outlet['e_PH'])

        # Case 2: No heat addition (Q <= 0)
        else:
            logging.warning("SolarTower: No solar heat input (Q <= 0).")
            self.E_P = np.nan
            self.E_F = np.nan

        #Case 3: Solar Tower with steam

        # Calculate exergetic efficiency
        self.epsilon = self.calc_epsilon()

        # Calculate exergy destruction 
        if not np.isnan(self.E_P):
            self.E_D = self.E_F - self.E_P
        else:
            self.E_D = self.E_F
