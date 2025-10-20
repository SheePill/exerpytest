import sys
import os
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Import the parser
from exerpy.parser.from_ebsilon.ebsilon_parser import EbsilonModelParser

# Define the path to the Ebsilon model file
model_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'csp.ebs'))

# Initialize the parser
parser = EbsilonModelParser(model_path, split_physical_exergy=True)
parser.initialize_model()
parser.simulate_model()

# Get a ThermoLiquid pipe
pipe = None
for j in range(1, parser.model.Objects.Count + 1):
    obj = parser.model.Objects.Item(j)
    if obj.IsKindOf(16):
        pipe_cast = parser.oc.CastToPipe(obj)
        if pipe_cast.FluidType == 20:
            pipe = pipe_cast
            break

print(f"Original pipe: {pipe.Name}")
print(f"  T={pipe.T.Value}°C, p={pipe.P.Value} bar")
print(f"  h={pipe.H.Value} kJ/kg, s={pipe.S.Value} kJ/kgK")

print("\n" + "="*80)
print("APPROACH: TEMPORARILY MODIFY PIPE TO GET h0, s0")
print("="*80)

# Save original values
T_original = pipe.T.Value
p_original = pipe.P.Value
h_original = pipe.H.Value
s_original = pipe.S.Value

# Set to Tmin
T_min = 200.0  # °C
p_ref = 1.01325  # bar (atmospheric pressure for reference state)

try:
    print(f"\nSetting pipe to T={T_min}°C, p={p_ref} bar")
    
    # Set new values
    pipe.T.Value = T_min
    pipe.P.Value = p_ref
    
    # Recalculate (might need to call calculate or similar)
    # Try to trigger recalculation
    try:
        parser.model.Calculate()
    except:
        pass
    
    # Read properties at reference state
    h0 = pipe.H.Value
    s0 = pipe.S.Value
    
    print(f"Reference state properties:")
    print(f"  h0 = {h0} kJ/kg")
    print(f"  s0 = {s0} kJ/kgK")
    
    # Calculate exergy
    T0_K = T_min + 273.15
    exergy_physical = (h_original - h0) - T0_K * (s_original - s0)
    
    print(f"\nPhysical exergy = (h - h0) - T0*(s - s0)")
    print(f"  = ({h_original} - {h0}) - {T0_K}*({s_original} - {s0})")
    print(f"  = {exergy_physical} kJ/kg")
    
finally:
    # Restore original values
    print(f"\nRestoring original values...")
    pipe.T.Value = T_original
    pipe.P.Value = p_original
    
    try:
        parser.model.Calculate()
    except:
        pass
    
    print(f"  T={pipe.T.Value}°C, p={pipe.P.Value} bar")
    print(f"  h={pipe.H.Value} kJ/kg, s={pipe.S.Value} kJ/kgK")
