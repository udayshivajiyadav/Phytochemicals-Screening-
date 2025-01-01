#takes the admet properties of all the molecules from IMPPAT and then generates a rank 



import pandas as pd

# Load the CSV file into a DataFrame
df = pd.read_csv('D:/Uday Paper/ADMET for IMPPAT Files/admet_properties.csv')

# Define Lipinski's Rule of Five ranges
mw_range = (0, 500)
logp_range = (0, 5)
tpsa_range = (0, 140)
h_donors_range = (0, 5)
h_acceptors_range = (0, 10)

# Check if each molecule meets each Lipinski criterion
df['Meets_MW'] = df['MolecularWeight'].between(*mw_range)
df['Meets_LogP'] = df['LogP'].between(*logp_range)
df['Meets_TPSA'] = df['TPSA'].between(*tpsa_range)
df['Meets_HDonors'] = df['NumHDonors'].between(*h_donors_range)
df['Meets_HAcceptors'] = df['NumHAcceptors'].between(*h_acceptors_range)

# Count the number of criteria met by each molecule
df['Criteria_Met'] = (df['Meets_MW'].astype(int) +
                      df['Meets_LogP'].astype(int) +
                      df['Meets_TPSA'].astype(int) +
                      df['Meets_HDonors'].astype(int) +
                      df['Meets_HAcceptors'].astype(int))

# Rank the molecules based on the number of criteria met
# 5 criteria met -> rank 1, 4 criteria met -> rank 2, and so on
df['Lipinski_Rank'] = df['Criteria_Met'].max() - df['Criteria_Met'] + 1

# Sort the DataFrame by Lipinski's Rank
df.sort_values(by='Lipinski_Rank', inplace=True)

# Add a new "rank" column
df.insert(0, 'rank', df['Lipinski_Rank'])

# Write the ranked DataFrame to a new CSV file
df.to_csv('ranked_molecules.csv', index=False)

print("Ranking saved to 'ranked_molecules.csv'")
