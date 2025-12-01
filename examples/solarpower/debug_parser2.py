import os

script_dir = os.path.dirname(__file__)
model_path = os.path.abspath(os.path.join(script_dir, 'Model2.ebs'))

from exerpy.parser.from_ebsilon.ebsilon_parser import EbsilonModelParser
import logging

# Enable detailed logging
logging.basicConfig(level=logging.INFO)

# Parse the model with detailed logging
parser = EbsilonModelParser(model_path, split_physical_exergy=False)
parser.initialize_model()
parser.simulate_model()
parser.parse_model()

print("\n\n===== DEBUG INFO =====")
print(f"Heatflux components to postprocess: {len(parser.heatflux_to_postprocess)}")
for hf in parser.heatflux_to_postprocess:
    print(f"  - Name: {hf['name']}, Q_Solar: {hf.get('Q_Solar')}, Unit: {hf.get('Q_Solar_unit')}")

print(f"\nAll components parsed:")
for group, comps in parser.components_data.items():
    print(f"  {group}:")
    for comp_name in comps.keys():
        print(f"    - {comp_name}")

print(f"\nConnections data keys: {list(parser.connections_data.keys())}")
