import os

script_dir = os.path.dirname(__file__)
model_path = os.path.abspath(os.path.join(script_dir, 'parabolheaders.ebs'))

from exerpy import ExergyAnalysis
from exerpy.components.nodes.collecting_header import CollectingHeader
from exerpy.components.nodes.distributing_header import DistributingHeader

ean = ExergyAnalysis.from_ebsilon(model_path, split_physical_exergy=False)


fuel = {"inputs": ['Oil'], "outputs": ['Oil_1']}
product = {"inputs": ['Water'], "outputs": ['Steam']}
loss = {"inputs": [], "outputs": []}



def print_header_details(analysis):
    """
    Print details of all collecting headers in the system.
    """
    for component in analysis.components.values():
        if isinstance(component, CollectingHeader):
            print(f"\nCollecting Header: {component.name}")
            print("Inlets:")
            for inlet_name, inlet in component.inl.items():
                print(f"  {inlet_name}:")
                print(f"    Mass flow: {inlet.get('m', 'N/A')} kg/s")
                print(f"    Temperature: {inlet.get('T', 'N/A')} K")
                print(f"    Pressure: {inlet.get('p', 'N/A')} bar")
                print(f"    Physical exergy: {inlet.get('e_PH', 'N/A')} kJ/kg")
            
            print("Outlets:")
            for outlet_name, outlet in component.outl.items():
                print(f"  {outlet_name}:")
                print(f"    Mass flow: {outlet.get('m', 'N/A')} kg/s")
                print(f"    Temperature: {outlet.get('T', 'N/A')} K")
                print(f"    Pressure: {outlet.get('p', 'N/A')} bar")
                print(f"    Physical exergy: {outlet.get('e_PH', 'N/A')} kJ/kg")
        elif isinstance(component, DistributingHeader):
            print(f"\nDistributing Header: {component.name}")
            print("Inlets:")
            for inlet_name, inlet in component.inl.items():
                print(f"  {inlet_name}:")
                print(f"    Mass flow: {inlet.get('m', 'N/A')} kg/s")
                print(f"    Temperature: {inlet.get('T', 'N/A')} K")
                print(f"    Pressure: {inlet.get('p', 'N/A')} bar")
                print(f"    Physical exergy: {inlet.get('e_PH', 'N/A')} kJ/kg")
            
            print("Outlets:")
            for outlet_name, outlet in component.outl.items():
                print(f"  {outlet_name}:")
                print(f"    Mass flow: {outlet.get('m', 'N/A')} kg/s")
                print(f"    Temperature: {outlet.get('T', 'N/A')} K")
                print(f"    Pressure: {outlet.get('p', 'N/A')} bar")
                print(f"    Physical exergy: {outlet.get('e_PH', 'N/A')} kJ/kg")

# Run the analysis
ean.analyse(E_F=fuel, E_P=product, E_L=loss)
ean.exergy_results()

# Print details of collecting headers
print("\n=== Collecting Header Details ===")
print_header_details(ean)