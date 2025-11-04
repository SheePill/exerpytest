import os

script_dir = os.path.dirname(__file__)
model_path = os.path.abspath(os.path.join(script_dir, '2.ebs'))

from exerpy import ExergyAnalysis

ean = ExergyAnalysis.from_ebsilon(model_path, split_physical_exergy=False)


fuel = {"inputs": ['Oil'], "outputs": ['Oil_1']}
product = {"inputs": ['Water'], "outputs": ['Steam']}
loss = {"inputs": [], "outputs": []}



ean.analyse(E_F=fuel, E_P=product, E_L=loss)
ean.exergy_results()