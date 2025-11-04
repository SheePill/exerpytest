import os

script_dir = os.path.dirname(__file__)
model_path = os.path.abspath(os.path.join(script_dir, 'Model3.ebs'))

from exerpy import ExergyAnalysis

ean = ExergyAnalysis.from_ebsilon(model_path, split_physical_exergy=False)

fuel = {"inputs": ['Oil_6'], "outputs": ['Oil_3']}
product = {"inputs": ['Electric'], "outputs": ['Electric_1', 'Electric_2', 'Electric_3']}
loss = {"inputs": ['Water_3'], "outputs": ['Water_4']}

ean.analyse(E_F=fuel, E_P=product, E_L=loss)
ean.exergy_results()