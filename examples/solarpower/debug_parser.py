import os
import json

script_dir = os.path.dirname(__file__)
model_path = os.path.abspath(os.path.join(script_dir, 'Model2.ebs'))

from exerpy.parser.from_ebsilon.ebsilon_parser import run_ebsilon

# Parse the model
parsed_data = run_ebsilon(model_path, split_physical_exergy=False)

# Print all connection names
print("All connections in the model:")
for conn_name in sorted(parsed_data['connections'].keys()):
    print(f"  - {conn_name}")

print("\n\nConnection details for 'Collector' (if exists):")
if 'Collector' in parsed_data['connections']:
    print(json.dumps(parsed_data['connections']['Collector'], indent=2))
else:
    print("  'Collector' not found in connections!")

print("\n\nHeatflux-like connections (containing 'heat' or 'solar'):")
for conn_name, conn_data in parsed_data['connections'].items():
    if 'heat' in str(conn_data).lower() or 'solar' in str(conn_data).lower() or conn_data.get('kind') == 'heat':
        print(f"  - {conn_name}: kind={conn_data.get('kind')}")
