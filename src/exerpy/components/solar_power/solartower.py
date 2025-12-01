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
        self.F = None

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
        if not hasattr(self, 'inl') or not hasattr(self, 'outl') or len(self.outl) != 1:
            msg = "SolarTower requires at least one inlet and exactly one outlet."
            logging.error(msg)
            raise ValueError(msg)
        
        if len(self.inl) < 1:
            msg = "SolarTower requires at least one inlet stream."
            logging.error(msg)
            raise ValueError(msg)

        # Extract inlet and outlet streams
        inlet = self.inl[0]
        outlet = self.outl[0]
        
        # Extract solar heat input from logic connection (inlet[1])
        if 1 in self.inl and self.inl[1] is not None:
            logic_inlet = self.inl[1]
            # Check if this is a heat connection with energy data
            if logic_inlet.get("kind") == "heat":
                # Try to get energy value from E or energy_flow fields
                self.F = logic_inlet.get("E") or logic_inlet.get("energy_flow")
                if self.F is None:
                    print(f"DEBUG: Warning - Logic connection missing energy data. Keys: {list(logic_inlet.keys())}")
        
        if self.F is None:
            print(f"DEBUG: Error - No valid solar heat input found, F is still None")

        # # Solar exergy calculations
        # # The sun's surface temperature [K]
        # T_SUN = 5778
        # #calculate Q reduction factor for exergy calculations
        # # E = Q * (1 - 4/3 * T_ambient/T_sun)
        # alpha = 1 - (4/3) * (T0 / T_SUN)
        # self.F = self.F
        
       # Case 1: Both inlet and outlet above ambient (normal operation)
        if split_physical_exergy:
                # Product is the increase in physical exergy
                self.E_P = outlet['m'] * (outlet['e_PH'] - inlet['e_PH'])
                self.E_F = abs(self.F)
        else:
                # Product is the increase in physical exergy
                self.E_P = outlet['m'] * (outlet['e_PH'] - inlet['e_PH'])
                self.E_F = abs(self.F)
                
        
        # Calculate exergetic efficiency
        self.epsilon = self.calc_epsilon()
        
        # Calculate exergy destruction 
        if not np.isnan(self.E_P):
            self.E_D = self.E_F - self.E_P
        else:
            self.E_D = self.E_F
        
        # Log the results
        logging.info(
            f"Solar Tower exergy balance calculated: "
            f"E_P={self.E_P:.2f}, E_F={self.E_F:.2f}, E_D={self.E_D:.2f}, "
            f"Efficiency={self.epsilon:.2%}"
        )    
