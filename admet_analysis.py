#carries out ADMET analysis of compounds

import csv
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

# Sample list of SMILES notations
df1 = pd.read_csv('D:/Uday Paper/ADMET for IMPPAT Files/admet_properties.csv')
smiles_list = df1['SMILES_form'].iloc[:17844].tolist()

def calculate_admet(smiles):
    # Convert SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        raise ValueError(f"Invalid SMILES notation: {smiles}")
    
    # Example properties calculation
    properties = {
        "SMILES": smiles,
        "MolecularWeight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),  # Topological Polar Surface Area
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol)
    }
    
    # Placeholder for actual ADMET prediction
    # You would typically use an external service here
    # For example, you could send a request to ADMETlab API
    
    admet_predictions = {
        "Absorption": "High",  # Placeholder value
        "Distribution": "Low",  # Placeholder value
        "Metabolism": "Low",  # Placeholder value
        "Excretion": "Medium",  # Placeholder value
        "Toxicity": "Low"  # Placeholder value
    }
    
    # Combine basic properties with ADMET predictions
    properties.update(admet_predictions)
    
    return properties

# Process each SMILES notation in the list
results = []
for smiles in smiles_list:
    try:
        properties = calculate_admet(smiles)
        results.append(properties)
    except ValueError as e:
        print(e)

# Define CSV file name
csv_file = 'admet_properties_test150.csv'

# Write the results to the CSV file
with open(csv_file, mode='w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)

print(f"Results have been written to {csv_file}")
