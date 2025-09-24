import logging

import numpy as np

from exerpy.components.component import Component
from exerpy.components.component import component_registry

@component_registry
class parabolic_trough(Component):

    def __init__(self, **kwargs):
        r"""
        Initializes the Heliostat component.
        """
        super().__init__(**kwargs)

    def calc_exergy_balance(self, T0: float, p0: float, split_physical_exergy) -> None:
        r"""
        Method to calculate the exergy balance of the Parabolic Trough component.
        """
        
        # Validate the number of inlets and outlets
        if not hasattr(self, 'inl') or not hasattr(self, 'outl') or len(self.inl) != 1 or len(self.outl) != 1:
            msg = "Parabolic Trough requires exactly one inlet and one outlet."
            logging.error(msg)
            raise ValueError(msg)
        
        # Extract inlet and outlet streams
        r"""
        component has no inlets or outlets?
        """