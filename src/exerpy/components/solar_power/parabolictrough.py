import logging

import numpy as np

from exerpy.components.component import Component
from exerpy.components.component import component_registry

@component_registry
class ParabolicTrough(Component):

    def __init__(self, **kwargs):
        r"""
        Initializes the Parabolic Trough component.

        Parameters
        ----------
        Q_Solar : float, optional
            Solar heat input from Ebsilon, [W]
        **kwargs : dict
            Additional keyword arguments passed to the parent class.
        """
        super().__init__(**kwargs)
        
        # Store the solar heat input from Ebsilon
        self.Q_Solar = kwargs.get('Q_Solar', None)
        if self.Q_Solar is None:
            logging.warning("Q_Solar not provided for Parabolic Trough component.")
        #self.Q_Loss = kwargs.get('QLOSS', None)
        #if self.Q_Loss is None:
            #logging.warning("QLOSS not provided for Parabolic Trough component.")
        self.Q_eff = kwargs.get('QEFF', None)
        if self.Q_eff is None:
            logging.warning("QEFF not provided for Parabolic Trough component.")
        
        # Extract NBRANCH value from heat flux connection and store as beta
        self.beta = kwargs.get('NBRANCH', None)
        if self.beta is None:
            logging.warning("NBRANCH not provided for Parabolic Trough component.")
        print(f"[DEBUG] Parabolic Trough '{self.name}' initialized with NBRANCH={self.beta}")
        

    def calc_exergy_balance(self, T0: float, p0: float, split_physical_exergy) -> None:
        r"""
        Method to calculate the exergy balance of the Parabolic Trough component.
        """
        
        # Validate the number of inlets and outlets
        if not hasattr(self, 'inl') or not hasattr(self, 'outl'):
            msg = "Collector requires two inlet connections and one outlet connection."
            logging.error(msg)
            raise ValueError(msg)
        
        if len(self.inl) != 1:
            msg = "Collector requires exactly two inlets."
            logging.error(msg)
            raise ValueError(msg)

        if len(self.outl) != 2:
            msg = "Collector requires exactly one outlet."
            logging.error(msg)
            raise ValueError(msg)
        

        # Extract inlet and outlet streams
        inlet = self.inl[0]
        outlet = self.outl[0]
        # Solar exergy calculations
        # The sun's surface temperature [K]
        T_SUN = 5778
        #calculate Q reduction factor for exergy calculations
        # E = Q * (1 - 4/3 * T_ambient/T_sun)
        alpha = 1 - (4/3) * (T0 / T_SUN)

        # Calculate exergy of incoming Heat
        self.E_Solar = self.Q_Solar * alpha * self.beta # Total incoming solar exergy
        #self.E_loss = self.Q_Loss * alpha    # Exergy losses from radiation
        self.E_eff = self.Q_eff * alpha   # Effective solar exergy

        # Case 1: Both inlet and outlet above ambient (normal operation)
        if inlet['T'] >= T0 and outlet['T'] >= T0:
            if split_physical_exergy:
                # Product is the increase in physical exergy
                self.E_P = self.beta * outlet['m'] * (outlet['e_PH'] - inlet['e_PH'])
                self.E_F = self.E_Solar
            else:
                # Product is the increase in physical exergy
                self.E_P = self.beta * outlet['m'] * (outlet['e_PH'] - inlet['e_PH'])
                #self.E_F = outlet['m'] * (outlet['e_PH'] - inlet['e_PH'])
                self.E_F = self.E_Solar
        

        # Calculate exergy destruction
        if np.isnan(self.E_P):
            self.E_D = self.E_F
        else:
            self.E_D = self.E_F - self.E_P 

        # Calculate exergy efficiency
        self.epsilon = self.calc_epsilon()

        # Log the results
        logging.info(
            f"Parabolic-Trough exergy balance calculated: "
            f"E_P={self.E_P:.2f}, E_F={self.E_F:.2f}, E_D={self.E_D:.2f}, "
            f"Efficiency={self.epsilon:.2%}"
        )
