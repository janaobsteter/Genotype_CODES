import pandas as pd
from scipy.stats.stats import pearsonr



scenarios = ["10_10", "50_50", "100_100", "0_100"]
rep = [4]

# Initiate the dataframe for TGVs for all bulls and bulls used in the home population
allBulls_final = pd.DataFrame()
homeFathersMean_final = pd.DataFrame()

# Check the results for each rep and each scenario within
for rep in rep:
    for scenario in scenarios:
        # Read in the pedigree file with the categories
        ped = pd.read_csv("GenGen" + str(rep) + "_" + scenario + "13/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
        # Read the file with individual and population (Group) information
        ps = pd.read_csv("PopulationSplit.txt", names = ['Group', 'Indiv'], header=1, low_memory=False)
        # Merge the files
        ped = ped.merge(ps, on="Indiv")
        # Compute the mean true genetic values for all genomically tested bulls in the home and import for trait 1 and trait 3
        allBulls = ped[ped.cat == "gpb"].groupby(["Group", "Generation"])[['gvNormUnres1', 'gvNormUnres3']].mean()
        allBulls.scenario = scenario
        allBulls.rep = rep
        # Concatenate the dataframes
        allBulls_final = allBulls_final.append(allBulls)


        # Extract only the fathers use in the home population and only genomically tested bulls (gen 40-60 in GenGen scenarios)
        homeFathers = ped[(ped.Indiv.isin(ped.Father[ped.Group == "home"])) & (ped.cat == "gpb")]
        # Compute the mean EBVs for import and home fathers used in the home population for trait 1 and trait 3
        homeFathersMean = homeFathers[homeFathers.cat == "gpb"].groupby(["Group", "Generation"])[['gvNormUnres1', 'gvNormUnres3']].mean()
        homeFathersMean.scenario = scenario
        homeFathersMean.rep = rep
        # Concatenate the dataframes
        homeFathersMean_final = homeFathersMean_final.append(homeFathersMean)



# Write the files
allBulls_final.to_csv("AllBulls_TGVs.csv")
# Write the file for the bulls used in home population
homeFathersMean_final.to_csv("HomeFathers_TGVs.csv")