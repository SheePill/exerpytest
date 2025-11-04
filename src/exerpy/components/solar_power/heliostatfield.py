import logging

import numpy as np

from exerpy.components.component import Component
from exerpy.components.component import component_registry

@component_registry
class Heliostatfield(Component):

    def __init__(self, **kwargs):
        r"""
        Initializes the Heliostat component.
        """
        super().__init__(**kwargs)
        self.Q_Solar = kwargs.get('Q_Solar', None)
        self.Q_eff = kwargs.get('RQINC', None)
    
    
    def calc_exergy_balance(self, T0: float, p0: float, split_physical_exergy) -> None:
        r"""
        Method to calculate the exergy balance of the Heliostat Field component.
        """
        # Validate the number of inlets and outlets
        if not hasattr(self, 'inl') or not hasattr(self, 'outl') or len(self.inl) != 1 or len(self.outl) != 1:
            msg = "Heliostat requires exactly one inlet and one outlet."
            logging.error(msg)
            raise ValueError(msg) 

        # Get the first inlet and outlet keys (they might not be 0-based)
        #inlet_key = list(self.inl.keys())[0]
        #outlet_key = list(self.outl.keys())[0]

        # Extract inlet and outlet streams
        #inlet = self.inl[inlet_key]
        #outlet = self.outl[outlet_key]

        # Solar exergy calculations
        # The sun's surface temperature [K]
        T_SUN = 5778
        #calculate Q reduction factor for exergy calculations
        # E = Q * (1 - 4/3 * T_ambient/T_sun)
        alpha = 1 - (4/3) * (T0 / T_SUN)

        # Calculate exergy of incoming Heat
        self.E_F = self.Q_Solar * alpha  # Total incoming solar exergy
        # Calculate exergy of effective Heat
        self.E_P = self.Q_eff * alpha  # Exergy of effective heat output
        #Calculate exergy destruction
        self.E_D = self.E_F - self.E_P

        # Assign exergy values to connections
        for logic_key in self.logic_outl:
            self.logic_outl[logic_key]['E'] = self.E_P  # Assign exergy to logic output
            
        # Calculate exergy efficiency
        self.epsilon = self.calc_epsilon()
        
        # Log the results
        logging.info(
            f"Heliostat Field Exergy Balance: "
            f"E_P={self.E_P:.2f}, E_F={self.E_F:.2f}, E_D={self.E_D:.2f}, "
            f"Efficiency={self.epsilon:.2%}"
        )