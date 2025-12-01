import os

script_dir = os.path.dirname(__file__)
model_path = os.path.abspath(os.path.join(script_dir, '4.5.ebs'))

from exerpy import ExergyAnalysis


ean = ExergyAnalysis.from_ebsilon(model_path, split_physical_exergy=False)


fuel = {"inputs": ['Heliostatenfeld_1_heatflux'], "outputs": []}
product = {"inputs": ['Electric'], "outputs": []}
loss = {"inputs": ['Water_2','Oil_3'], "outputs": ['Water_1','Oil']}


ean.analyse(E_F=fuel, E_P=product, E_L=loss)
ean.exergy_results()