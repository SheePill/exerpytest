import logging
from typing import Any
from typing import Optional

from CoolProp.CoolProp import PropsSI as CP

from . import __ebsilon_available__
from .utils import EpGasTableStub
from .utils import EpSteamTableStub
from .utils import require_ebsilon

# Import Ebsilon classes if available
if __ebsilon_available__:
    from EbsOpen import EpGasTable
    from EbsOpen import EpSteamTable
else:
    EpSteamTable = EpSteamTableStub
    EpGasTable = EpGasTableStub

from exerpy.functions import convert_to_SI

from .ebsilon_config import substance_mapping
from .ebsilon_config import unit_id_to_string


@require_ebsilon
def calc_X_from_PT(app: Any, pipe: Any, property: str, pressure: float, temperature: float) -> Optional[float]:
    """
    Calculate a thermodynamic property (enthalpy or entropy) for a given stream based on pressure and temperature.

    This method takes pressure and temperature values and calculates the specified property for any fluid stream.
    It automatically handles the composition of the stream by setting up appropriate fluid properties and analysis
    parameters based on the stream's fluid type and composition.

    Parameters
    ----------
    app : Ebsilon application instance
        The Ebsilon application used for creating fluid and analysis objects.
    pipe : Stream object
        The stream object containing fluid and composition information.
    property : str
        The thermodynamic property to calculate, either 'H' for enthalpy or 'S' for entropy.
    pressure : float
        The pressure value (in bar).
    temperature : float
        The temperature value (in °C).

    Returns
    -------
    float
        The calculated value of the specified property (in J/kg for enthalpy, J/kgK for entropy).
        Returns None if an invalid property is specified or an error occurs during calculation.

    Raises
    ------
    ValueError
        If Ebsilon returns the sentinel value -999.0 indicating a failed calculation.
    Exception
        Logs an error and returns None if any other exception occurs during property calculation.
    """

    # Create a new FluidData object
    fd = app.NewFluidData()

    # Retrieve the fluid type from the stream
    fd.FluidType = (pipe.Kind-1000)

    if fd.FluidType == 3 or fd.FluidType == 4:  # steam or water
            t_sat = CP('T', 'P', pressure, 'Q', 0, 'water')
            if temperature > t_sat:
                fd.FluidType = 3  # steam
                fd.SteamTable = EpSteamTable.epSteamTableFromSuperiorModel
                fdAnalysis = app.NewFluidAnalysis()
            else:
                fd.FluidType == 4  # water
                fdAnalysis = app.NewFluidAnalysis()

    elif fd.FluidType == 15:  # 2PhaseLiquid
        fd.Medium = pipe.FMED.Value
        fdAnalysis = app.NewFluidAnalysis()

    elif fd.FluidType == 16:  # 2PhaseGaseous
        fd.Medium = pipe.FMED.Value
        fdAnalysis = app.NewFluidAnalysis()

    elif fd.FluidType == 17:  # Salt water
        fd.Medium = pipe.FMED.Value
        fdAnalysis = app.NewFluidAnalysis()

    elif fd.FluidType == 20:  # ThermoLiquid
        fd.Medium = pipe.FMED.Value
        fdAnalysis = app.NewFluidAnalysis()
        
    else:  # flue gas, air etc.
        fd.GasTable = EpGasTable.epGasTableFromSuperiorModel

        # Set up the fluid analysis based on stream composition
        fdAnalysis = app.NewFluidAnalysis()

        # Iterate through the substance_mapping and get the corresponding value from the pipe
        for substance_key, ep_substance_id in substance_mapping.items():
            fraction = getattr(pipe, substance_key).Value  # Dynamically access the fraction
            if fraction > 0:  # Only set substances with non-zero fractions
                fdAnalysis.SetSubstance(ep_substance_id, fraction)

    # Set the analysis in the FluidData object
    fd.SetAnalysis(fdAnalysis)

    # Validate property input
    if property not in ['S', 'H']:
        logging.error('Invalid property selected. You can choose between "H" (enthalpy) and "S" (entropy).')
        return None

    try:
        # Calculate the property based on the input property type
        if property == 'S':  # Entropy
            res = fd.PropertyS_OF_PT(pressure * 1e-5, temperature - 273.15)  # Ebsilon works with °C and bar
            res_SI = res * 1e3  # Convert kJ/kgK to J/kgK
        elif property == 'H':  # Enthalpy
            res = fd.PropertyH_OF_PT(pressure * 1e-5, temperature - 273.15)  # Ebsilon works with °C and bar
            res_SI = res * 1e3  # Convert kJ/kg to J/kg
        
        if res == -999.0:
            raise ValueError(
                f"Calculation with Ebsilon property failed: {property} = -999.0 "
                f"for fluid {fd.Medium} at p={pressure} Pa and T={temperature} K. "
                f"It may helpful to set split_physical_exergy=False in the ExergyAnalysis constructor."
            )

        return res_SI

    except Exception as e:
        if isinstance(e, ValueError):
            raise
        logging.error(f"An error occurred during property calculation: {e}")
        return None


@require_ebsilon
def calc_eT(app: Any, pipe: Any, pressure: float, Tamb: float, pamb: float) -> float:
    """
    Calculate the thermal component of physical exergy.

    Parameters
    ----------
    app : Ebsilon application instance
        The Ebsilon application instance.
    pipe : Stream object
        The stream object containing thermodynamic properties.
    pressure : float
        The pressure value (in bar).
    Tamb : float
        The ambient temperature (in K).
    pamb : float
        The ambient pressure (in Pa).

    Returns
    -------
    float
        The thermal exergy component (in J/kg).
    """
    h_i = convert_to_SI('h', pipe.H.Value, unit_id_to_string.get(pipe.H.Dimension, "Unknown"))  # in SI unit [J / kg]
    s_i = convert_to_SI('s', pipe.S.Value, unit_id_to_string.get(pipe.S.Dimension, "Unknown"))  # in SI unit [J / kgK]
    h_A = calc_X_from_PT(app, pipe, 'H', pressure, Tamb)  # in SI unit [J / kg]
    s_A = calc_X_from_PT(app, pipe, 'S', pressure, Tamb)  # in SI unit [J / kgK]
    eT = h_i - h_A - Tamb * (s_i - s_A)  # in SI unit [J / kg]

    return eT


@require_ebsilon
def calc_eM(app: Any, pipe: Any, pressure: float, Tamb: float, pamb: float) -> float:
    """
    Calculate the mechanical component of physical exergy.

    Parameters
    ----------
    app : Ebsilon application instance
        The Ebsilon application instance.
    pipe : Stream object
        The stream object containing thermodynamic properties.
    pressure : float
        The pressure value (in bar).
    Tamb : float
        The ambient temperature (in K).
    pamb : float
        The ambient pressure (in Pa).

    Returns
    -------
    float
        The mechanical exergy component (in J/kg).
    """
    eM = convert_to_SI('e', pipe.E.Value, unit_id_to_string.get(pipe.E.Dimension, "Unknown")) - calc_eT(app, pipe, pressure, Tamb, pamb)

    return eM


def calc_eph_from_min(app: Any, pipe: Any, pressure: float, temperature: float) -> Optional[float]:
    """
    Calculate physical exergy using the minimum-valid-temperature reference state.

    Formula:
        E_PH = H - H_min - T_min * (S - S_min)

    Parameters
    ----------
    app : Ebsilon application instance
    pipe : Stream object (pipe_cast)
    pressure : float
        Stream pressure in SI [Pa]
    temperature : float
        Stream temperature in SI [K]

    Returns
    -------
    float or None
        Physical exergy [J/kg], or None if it cannot be computed.

    Notes
    -----
    - Supports both regular FluidData and UniversalFluidData (for UniversalSubstances like molten salts)
    - Automatically detects fluid type via pipe.HasUniversalFluidData property
    - If PropertyTMIN_OF_P returns -999.0, fallback search finds minimum valid temperature
    """
    import logging

    # 1) Determine if pipe uses UniversalFluidData or regular FluidData
    try:
        if hasattr(pipe, 'HasUniversalFluidData') and pipe.HasUniversalFluidData:
            # Use the UniversalFluidData already configured in the pipe
            fd = pipe.UniversalFluidData()
            logging.info(f"Using UniversalFluidData for {getattr(pipe, 'Name', '<unknown>')}")
            
        else:
            # Build regular FluidData as before
            fluid_type = (pipe.Kind - 1000)
            fd = app.NewFluidData()
            fd.FluidType = fluid_type
            fdAnalysis = app.NewFluidAnalysis()

            # Steam/Water handling
            if fluid_type == 3 or fluid_type == 4:  # steam or water
                try:
                    t_sat = CP('T', 'P', pressure, 'Q', 0, 'water')
                except Exception:
                    t_sat = None

                if t_sat is not None and temperature > t_sat:
                    fd.FluidType = 3  # steam
                    fd.SteamTable = EpSteamTable.epSteamTableFromSuperiorModel
                else:
                    fd.FluidType = 4  # water

            elif fluid_type in (15, 16, 17, 20):  # 2PhaseLiquid, 2PhaseGaseous, Salt water, ThermoLiquid
                if hasattr(pipe, 'FMED'):
                    fd.Medium = pipe.FMED.Value

            else:
                # flue gas, air etc.
                fd.GasTable = EpGasTable.epGasTableFromSuperiorModel
                # Set composition from mapping
                for substance_key, ep_substance_id in substance_mapping.items():
                    if hasattr(pipe, substance_key):
                        fraction = getattr(pipe, substance_key).Value
                        if fraction and fraction > 0:
                            fdAnalysis.SetSubstance(ep_substance_id, fraction)

            fd.SetAnalysis(fdAnalysis)
            
    except Exception as e:
        logging.error(f"Failed to prepare FluidData for {getattr(pipe, 'Name', '<unknown>')}: {e}")
        return None

    # 2) Read current state from the pipe (convert to SI)
    try:
        h_i = convert_to_SI('h', pipe.H.Value, unit_id_to_string.get(pipe.H.Dimension, "Unknown"))  # J/kg
        s_i = convert_to_SI('s', pipe.S.Value, unit_id_to_string.get(pipe.S.Dimension, "Unknown"))  # J/kgK
    except Exception as e:
        logging.error(f"Could not read H/S from pipe {getattr(pipe, 'Name', '<unknown>')}: {e}")
        return None

    # 3) Units for property calls
    p_bar = pressure * 1e-5       # Pa -> bar (Ebsilon expects bar)
    T_C   = temperature - 273.15  # K  -> °C

    # 4) Try the direct property first
    t_min_C = None
    try:
        tmin_candidate = fd.PropertyTMIN_OF_P(p_bar)  # returns °C or -999.0
        if tmin_candidate != -999.0:
            t_min_C = tmin_candidate
    except Exception as e:
        logging.debug(f"PropertyTMIN_OF_P raised: {e}")

    # 5) Fallback: probe to find the lowest valid temperature
    def is_valid_at(temp_C: float) -> bool:
        """Check if H(T,p) is valid at this temperature (not -999.0)."""
        h_kJ = fd.PropertyH_OF_PT(p_bar, temp_C)
        return h_kJ != -999.0

    if t_min_C is None:
        # Start near the actual stream temperature if valid; otherwise try to climb up to a valid point
        start_C = T_C if is_valid_at(T_C) else None
        if start_C is None:
            # Climb upward until we find a valid point (handles cases where T is below range)
            probe = T_C
            for _ in range(60):  # up to ~600°C span with 10°C steps
                probe += 10.0
                if is_valid_at(probe):
                    start_C = probe
                    break
        if start_C is None:
            logging.error(f"Could not find any valid temperature at p={p_bar:.4f} bar for {getattr(pipe, 'Name', '<unknown>')}")
            return None

        # Step downward until we hit invalid -> bracket found
        step = 20.0
        T_valid = start_C
        T_probe = start_C
        while is_valid_at(T_probe) and T_probe > -150.0:
            T_valid = T_probe
            T_probe -= step

        if is_valid_at(T_probe):
            # Never became invalid within scan range
            logging.warning(f"Did not find lower invalid bound; using lowest scanned valid {T_valid:.3f} °C as t_min for {getattr(pipe, 'Name', '<unknown>')}")
            t_min_C = T_valid
        else:
            T_invalid = T_probe
            # Binary refine between [T_invalid, T_valid] to the boundary
            tol = 0.05  # °C tolerance
            while (T_valid - T_invalid) > tol:
                mid = 0.5 * (T_valid + T_invalid)
                if is_valid_at(mid):
                    T_valid = mid
                else:
                    T_invalid = mid
            t_min_C = T_valid  # lowest valid temperature

    # 6) Compute h_min and s_min at T_min (slightly inside the valid region)
    epsilon = 1e-3
    tries = 0
    h_min_kJ = s_min_kJ = -999.0
    while tries < 5 and (h_min_kJ == -999.0 or s_min_kJ == -999.0):
        T_query = t_min_C + tries * max(epsilon, 0.05)
        h_min_kJ = fd.PropertyH_OF_PT(p_bar, T_query)
        s_min_kJ = fd.PropertyS_OF_PT(p_bar, T_query)
        tries += 1

    if h_min_kJ == -999.0 or s_min_kJ == -999.0:
        logging.error(
            f"Failed to evaluate h_min/s_min around T_min≈{t_min_C:.3f} °C "
            f"for {getattr(pipe, 'Name', '<unknown>')} at p={p_bar:.4f} bar."
        )
        return None

    # 7) Assemble the exergy
    h_min = h_min_kJ * 1e3  # kJ/kg -> J/kg
    s_min = s_min_kJ * 1e3  # kJ/kgK -> J/kgK
    t_min_K = (t_min_C + 273.15)

    e_ph = h_i - h_min - t_min_K * (s_i - s_min)

    logging.info(
        f"e_PH(min) for {getattr(pipe, 'Name', '<unknown>')}: "
        f"T_min={t_min_C:.2f} °C, h_min={h_min:.1f} J/kg, s_min={s_min:.3f} J/kgK, e_PH={e_ph:.1f} J/kg"
    )
    return e_ph



@require_ebsilon
def debug_substance_detection(app, pipe):
    """
    Debug helper to understand what substance/medium Ebsilon is using for a pipe.
    
    Parameters
    ----------
    app : IApplication
        Ebsilon application instance
    pipe : IPipe
        Pipe object to debug
        
    Returns
    -------
    dict
        Information about the substance detection
    """
    fluid_type_index = {
        1: "IAPWS-IF97 (Water/Steam)",
        3: "Ideal gas",
        4: "Real gas",
        6: "Numeric value",
        9: "Binary mixture",
        10: "Power/Mechanical",
        13: "Logic",
        20: "ThermoLiquid (substance 120)",
        1122: "Molten Salt (substance 122)"
    }
    
    result = {
        'pipe_name': pipe.Name,
        'fluid_type': pipe.FluidType,
        'fluid_type_name': fluid_type_index.get(pipe.FluidType, "Unknown"),
        'FMED': pipe.FMED.Value,
        'T': pipe.T.Value,
        'p': pipe.P.Value,
        'substance_attempts': {}
    }
    
    # Try to get substance properties
    try:
        # Try substance 120 (default for ThermoLiquid)
        sub120 = app.GetSubstance(120)
        result['substance_attempts']['120'] = {
            'exists': True,
            'name': sub120.Name if hasattr(sub120, 'Name') else 'N/A'
        }
    except Exception as e:
        result['substance_attempts']['120'] = {
            'exists': False,
            'error': str(e)
        }
    
    try:
        # Try substance 122 (Molten Salt)
        sub122 = app.GetSubstance(122)
        result['substance_attempts']['122'] = {
            'exists': True,
            'name': sub122.Name if hasattr(sub122, 'Name') else 'N/A'
        }
    except Exception as e:
        result['substance_attempts']['122'] = {
            'exists': False,
            'error': str(e)
        }
    
    # Try to change FMED and see what happens
    try:
        original_fmed = pipe.FMED.Value
        pipe.FMED.Value = 122  # Try to set to molten salt
        result['fmed_change_test'] = {
            'original': original_fmed,
            'attempted': 122,
            'success': pipe.FMED.Value == 122,
            'actual_value': pipe.FMED.Value
        }
        # Restore original
        pipe.FMED.Value = original_fmed
    except Exception as e:
        result['fmed_change_test'] = {
            'error': str(e)
        }
    
    return result


def _get_stream_properties(pipe, T0, p0):
    """
    Get thermodynamic properties from an Ebsilon pipe.
    For molten salt and other ThermoLiquid, use Ebsilon's calculated values.
    """
    properties = {}
    
    # Get current state from pipe
    properties['T'] = pipe.T.Value  # °C
    properties['p'] = pipe.P.Value  # bar
    properties['m'] = pipe.M.Value  # kg/s
    
    # For ThermoLiquid (molten salt, thermoil, etc.), Ebsilon calculates h and s
    if pipe.FluidType == 20:
        # Use Ebsilon's enthalpy and entropy directly
        properties['h'] = pipe.H.Value  # kJ/kg
        properties['s'] = pipe.S.Value  # kJ/kgK
        
        # Reference state properties - need to be calculated at T0, p0
        # We'll handle this separately
        
    elif pipe.FluidType == 1:  # Water/Steam
        properties['h'] = pipe.H.Value
        properties['s'] = pipe.S.Value
    
    return properties
