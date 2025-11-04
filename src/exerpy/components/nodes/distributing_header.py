import logging

import numpy as np

from exerpy.components.component import Component
from exerpy.components.component import component_registry

@component_registry
class DistributingHeader(Component):

    def __init__(self, **kwargs):
        r"""
        Initializes the Distributing Header component.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed to the parent class.
        """
        super().__init__(**kwargs)
        #Gets the number of branches if provided in Ebsilon model
        self.NBRANCH = kwargs.get('NBRANCH', None)
        if self.NBRANCH is None:
            logging.warning("NBRANCH not provided for Distributing Header component.")

    def calc_exergy_balance(self, T0: float, p0: float, split_physical_exergy) -> None:
        r"""
        Method to calculate the exergy balance of the Distributing Header component.
        """
        
        # Validate the number of inlets and outlets
        if not hasattr(self, 'inl') or not hasattr(self, 'outl') or len(self.inl) != 1 or len(self.outl) < 1:
            msg = "Distributing Header requires exactly one inlet and at least one outlet."
            logging.error(msg)
            raise ValueError(msg)
        

        outlet_list = list(self.outl.values())
        inlet_list = list(self.inl.values())
        
        # Setting x to NBRANCH attribute if provided
        # Adjust mass flow and heat flow of outlets based on NBRANCH so that total flow is conserved
        x = self.NBRANCH
        if len(self.outl) == 2 and x != 1.0:
            for outlet in self.outl.values():
                if 'm' in outlet:
                    outlet['m'] = outlet.get('m', 0) * x
                if 'Q' in outlet:
                    outlet['Q'] = outlet.get('Q', 0) * x

        
        E_in = sum(inlet.get("m", 0) * inlet.get("e_PH") for inlet in inlet_list)
        E_out = sum(outlet.get("m", 0) * outlet.get("e_PH") for outlet in outlet_list)
        self.E_P = np.nan
        self.E_F = np.nan
        self.E_D = E_in - E_out
        self.epsilon = np.nan

        # Log the results.
        logging.info(
            f"Exergy balance of Collecting-Header {self.name} calculated: "
            f"E_P={self.E_P:.2f}, E_F={self.E_F:.2f}, E_D={self.E_D:.2f}, "
            f"Efficiency={self.epsilon:.2%}"
        )