import logging

import numpy as np

from exerpy.components.component import Component
from exerpy.components.component import component_registry
from exerpy.functions import fluid_property_data

#anstatt Variable Connection heatflux benutzen damit es eleganter ist

@component_registry
class Heliostatfield(Component):

    def __init__(self, **kwargs):
        r"""
        Initializes the Heliostat component.
        """
        super().__init__(**kwargs)
        self.Q_Solar = kwargs.get('Q_Solar', None)
        self.F = None
        self.P = None
    
    
    def calc_exergy_balance(self, T0: float, p0: float, split_physical_exergy) -> None:
        r"""
        Method to calculate the exergy balance of the Heliostat Field component.
        """
        # Validate the number of inlets and outlets
        if not hasattr(self, 'outl') or not hasattr(self, 'outl'):
            msg = "Heliostat requires inlet and outlet connections."
            logging.error(msg)
            raise ValueError(msg)
        
        if len(self.outl) > 2: #FIXME: Should this be exactly one inlet?
            msg = "Heliostat requires exactly one inlet."
            logging.error(msg)
            raise ValueError(msg)
            
        if len(self.outl) > 2: #FIXME: Should this be exactly two outlets?
            msg = "Heliostat requires exactly two outlets."
            logging.error(msg)
            raise ValueError(msg)
        
        # Check if there are any outlets at all
        if len(self.outl) == 0:
            logging.warning(f"Heliostat Field '{self.name}' has no outlet connections.")
            self.E_F = 0
            self.E_P = 0
            self.E_D = 0
            self.epsilon = np.nan
            return
        
        # Extract inlet and outlet streams
        # Normalize placeholder zeros to None for safety (some parsers use 0 as placeholder)
        if hasattr(self, 'inl') and 0 in self.inl and self.inl[0] == 0:
            self.inl[0] = None
        if hasattr(self, 'outl') and 0 in self.outl and self.outl[0] == 0:
            self.outl[0] = None

        if 0 in self.outl and self.outl[0] is not None:
            logic_outlet = self.outl[0]
            # Check if this is a heat connection with energy data
            if logic_outlet.get("kind") == "heat":
                # Try to get energy value from E or energy_flow fields
                self.P = logic_outlet.get("E") or logic_outlet.get("energy_flow")
                if self.P is None:
                    print(f"DEBUG: Warning - Logic connection missing energy data. Keys: {list(logic_outlet.keys())}")
        

        
        # Solar exergy calculations     
        #calculate Q reduction factor for exergy calculations
        # E = Q * (1 - 4/3 * T_ambient/T_sun)
        # The sun's surface temperature [K]
        T_SUN = 5778
        alpha = 1 - (4/3) * (T0 / T_SUN)

        # Calculate exergy of incoming Heat
        self.E_F = self.Q_Solar * alpha  # Total incoming solar exergy
        # Calculate exergy of effective Heat
        # Use 'energy_flow' key for heat connections from the parser
        self.E_P = self.P * alpha # Exergy of effective heat output
        #Calculate exergy destruction
        self.E_D = self.E_F - self.E_P

        # Update the connection data dictionary (if present) so the rest of the pipeline
        try:
            if logic_outlet is not None:
                # store exergy flow in W and set unit to heat SI unit
                logic_outlet["E"] = self.E_P
                logic_outlet["E_unit"] = fluid_property_data["heat"]["SI_unit"]
        except Exception:
            # defensive: don't let connection update failures break the exergy calculation
            logging.debug(f"Could not write back exergy to connection for Heliostat '{self.name}'")
        try:
            if logic_inlet is not None:
                # store exergy flow in W and set unit to heat SI unit
                logic_inlet["E"] = self.E_F
                logic_inlet["E_unit"] = fluid_property_data["heat"]["SI_unit"]
        except Exception:
            # defensive: don't let connection update failures break the exergy calculation
            logging.debug(f"Could not write back exergy to connection for Heliostat '{self.name}'")


        # Calculate exergy efficiency
        self.epsilon = self.calc_epsilon()
        
        # Log the results
        logging.info(
            f"Heliostat Field Exergy Balance: "
            f"E_P={self.E_P:.2f}, E_F={self.E_F:.2f}, E_D={self.E_D:.2f}, "
            f"Efficiency={self.epsilon:.2%}"
        )