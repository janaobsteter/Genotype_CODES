#!/usr/bin/env python
# coding: utf-8


import cyvcf2
import tsinfer
import os
import pandas as pd
import json
import sys
#import tqdm

vcfFile = sys.argv[1]
vcfFileName = vcfFile.split("/")[-1]
sampleFile = str(vcfFileName) + ".samples"
#contigsize = sys.argv[2]

# Add progress bar
#progressbar = tqdm.tqdm(total=samples.sequence_length, desc="Read VCF", unit='bp')

# Define the functions to reading from a vcf
def add_diploid_sites(vcf, samples):
    """
    Read the sites in the vcf and add them to the samples object, reordering the
    alleles to put the ancestral allele first, if it is available.
    """
        
    pos = 0
    for variant in vcf:  # Loop over variants, each assumed at a unique site
        # Update on the progress
        progressbar.update(variant.POS - pos)
        
        if pos == variant.POS:
            raise ValueError("Duplicate positions for variant at position", pos)
        else:
            pos = variant.POS
        if any([not phased for _, _, phased in variant.genotypes]):
            raise ValueError("Unphased genotypes for variant at position", pos)
        alleles = [variant.REF] + variant.ALT
        ancestral = variant.INFO.get("AA", variant.REF)
        # Ancestral state must be first in the allele list.
        ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
        allele_index = {
            old_index: ordered_alleles.index(allele)
            for old_index, allele in enumerate(alleles)
        }
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = [
            allele_index[old_index]
            for row in variant.genotypes
            for old_index in row[0:2]
        ]
        samples.add_site(pos, genotypes=genotypes, alleles=alleles)
        
def chromosome_length(vcf):
    assert len(vcf.seqlens) == 1
    return vcf.seqlens[0]

def add_diploid_individuals(vcf, samples, populations):
    for name, population in zip(vcf.samples, populations):
        samples.add_individual(ploidy=2, metadata={"name": name}, population=population)



# Read in the vcf file
vcf = cyvcf2.VCF(vcfFile)

# Read in the meta information
meta = pd.read_csv("Meta.csv")
# Retain the information about the ID and the subspecie
meta = meta[['ID', 'Type']]
metaPop = [list(meta.Type[meta.ID == x])[0] for x in vcf.samples]

# This is a modified add populations function for drones
def add_populations_drones(vcf, samples, metaData):
    """
    Add tsinfer Population objects for drone subspecies, stored in metaPop object, and return a list od IDs correposning to the VCF samples
    Input: vcf = cyvcf2 VCF object, samples is tsinfer SampleData and metaData is a pandas dataframe with id in ID column and populations in Type column
    Return: A list of population indexes
    """
    pop_lookup = {}
    pop_lookup["carnica"] = samples.add_population(metadata={"subspecies": "carnica", "lineage": "C"})
    pop_lookup["mellifera"] = samples.add_population(metadata={"subspecies": "mellifera", "lineage": "M"})
    pop_lookup["scutellata"] = samples.add_population(metadata={"subspecies": "scutellata", "lineage": "A"})
    pop_lookup["caucasica"] = samples.add_population(metadata={"subspecies": "caucasica", "lineage": "O"})
    pop_lookup["unicolor"] = samples.add_population(metadata={"subspecies": "unicolor", "lineage": "A"})
    pop_lookup["hybrid"] = samples.add_population(metadata={"subspecies": "hybrid", "lineage": "H"})
    pop_lookup["capensis"] = samples.add_population(metadata={"subspecies": "capensis", "lineage": "A"})
    pop_lookup["ligustica"] = samples.add_population(metadata={"subspecies": "ligustica", "lineage": "C"})
    return [pop_lookup[x] for x in [list(meta.Type[meta.ID == x])[0] for x in vcf.samples]]


# Create the samples
print("Do the inference")
with tsinfer.SampleData(
    path="Drones.samples", sequence_length=250e06
) as samples:
    populations = add_populations_drones(vcf, samples, meta)
    add_diploid_individuals(vcf, samples, populations)
    add_diploid_sites(vcf, samples)

# Print the summary of the samples created
print(
    "Sample file created for {} samples ".format(samples.num_samples)
    + "({} individuals) ".format(samples.num_individuals)
    + "with {} variable sites.".format(samples.num_sites),
    flush=True,
)

# Do the inference on the 10 SNPs
drone_ts = tsinfer.infer(samples)
print(
    "Inferred tree sequence `{}`: {} trees over {} Mb".format(
        "drone_ts", drone_ts.num_trees, drone_ts.sequence_length / 1e6
    )
)


# Check the metadata
for sample_node_id in drone_ts.samples():
    individual_id = drone_ts.node(sample_node_id).individual
    population_id = drone_ts.node(sample_node_id).population
    print(
        "Node",
        sample_node_id,
        "labels genome sampled from",
        json.loads(drone_ts.individual(individual_id).metadata),
        "in",
        json.loads(drone_ts.population(population_id).metadata)["subspecies"],
    )

    print("Draw a tree")
# Set the colours of the subspecies for plotting    
colours = {"hybrid": "grey", "carnica": "yellow", "ligustica": "orange", "mellifera": "red", "caucasica": "blue", "scutellata": "black", "unicolor": "blue", "capensis": "purple", }
# Set the colours for the individuals
colours_for_node = {}

# If you want to plot by lineages
# Set the colours of the lineages for plotting
coloursL = {"A": "black", "C": "yellow", "M": "red", "O": "blue", "H": "grey"}
colours_for_nodeL = {}

for n in drone_ts.samples():
    population_data = drone_ts.population(drone_ts.node(n).population)
    colours_for_node[n] = colours[json.loads(population_data.metadata)["subspecies"]]
    colours_for_nodeL[n] = coloursL[json.loads(population_data.metadata)["lineage"]]

individual_for_node = {}
for n in drone_ts.samples():
    individual_data = drone_ts.individual(drone_ts.node(n).individual)
    individual_for_node[n] = json.loads(individual_data.metadata)["name"]
    
    
tree = drone_ts.at(5e6)
tree.draw(
    path="tree_at_1Mb.svg",
    height=700,
    width=1200,
    node_labels=individual_for_node,
    node_colours=colours_for_node
)

tree.draw(
    path="tree_at_1Mb_lineages.svg",
    height=700,
    width=1200,
    node_labels=individual_for_node,
    node_colours=colours_for_nodeL
)


# In[ ]:


#Obtain node info - but just for samples!
# There is also intermediate nodes!!! (total #nodes > #nodes for samples)
sample_nodes = [drone_ts.node(n) for n in drone_ts.samples()]

#Get samples ids
sample_ids = [n.id for n in sample_nodes]

# Get sample names
sample_names = [
    json.loads(drone_ts.individual(n.individual).metadata)["name"]
    for n in sample_nodes
]
# Get sample population
sample_pops = [
    json.loads(drone_ts.population(n.population).metadata)["lineage"]
    for n in sample_nodes
]
# #print(sample_pops)
# print(len(sample_pops))

# This one takes a long time since there is many nodes!!!
# node_names = [
#     json.loads(drone_ts.individual(n.individual).metadata)["name"]
#     if n in sample_nodes else None for n in drone_ts.nodes()
# ]

# Create a dictionary holding the sample node ids within each population (in this case sample_pops are lineages)
pops = dict()
for (node, pop) in zip(sample_ids, sample_pops):
    if not pop in pops.keys():
        pops[pop] = [node]
    else:
        pops[pop].append(node)



# In[ ]:


print("Do the statistics.")
# Obtain Fst and write a table   
Fst = drone_ts.Fst([pops["A"], pops["C"], pops["M"], pops["O"]], indexes=[(0,1), (0,2), (0, 3), (1, 2), (1, 3), (2,3)])
#print(Fst)
fstDF = pd.DataFrame(Fst, columns = ["Fst"])
fstDF.loc[:, "Group1"] = ["A", "A", "A", "C", "C", "M"]
fstDF.loc[:, "Group2"] = ["C", "M", "O", "M", "O", "O"]
fstDF.to_csv("FST.csv", index=None)


# In[ ]:


# Obtain Tajima's D and write a table   
tajimaD = drone_ts.Tajimas_D([pops["A"], pops["C"], pops["M"], pops["O"], pops["H"]])
#print(tajimaD)
tajimaDF = pd.DataFrame(tajimaD, columns = ["TajimasD"])
tajimaDF.loc[:, "Group"] = ["A", "C", "M", "O", "H"]
tajimaDF.to_csv("TajimasD.csv", index=None)


# In[ ]:


# Divergence 
divergence = drone_ts.divergence([pops["A"], pops["C"], pops["M"], pops["O"]], indexes=[(0,1), (0,2), (0, 3), (1, 2), (1, 3), (2,3)])
#print(divergence)
divergenceDF = pd.DataFrame(divergence, columns = ["Divergence"])
divergenceDF.loc[:, "Group1"] = ["A", "A", "A", "C", "C", "M"]
divergenceDF.loc[:, "Group2"] = ["C", "M", "O", "M", "O", "O"]
divergenceDF.to_csv("Divergence.csv", index=None)


# In[ ]:


# Diversity and write a table   
diversity = drone_ts.diversity([pops["A"], pops["C"], pops["M"], pops["O"], pops["H"]])
diversityDF = pd.DataFrame(diversity, columns = ["TajimasD"])
diversityDF.loc[:, "Group"] = ["A", "C", "M", "O", "H"]
diversityDF.to_csv("Diversity.csv", index=None)


# In[ ]:


# Allele frequency spectrum 
np.savetxt("AFS_Alineage.csv", drone_ts.allele_frequency_spectrum([pops["A"]]))
np.savetxt("AFS_Mlineage.csv", drone_ts.allele_frequency_spectrum([pops["M"]]))
np.savetxt("AFS_Clineage.csv", drone_ts.allele_frequency_spectrum([pops["C"]]))
np.savetxt("AFS_Olineage.csv", drone_ts.allele_frequency_spectrum([pops["O"]]))


# In[ ]:


# f2, f3, f4 statistics 
f2 = drone_ts.f2([pops["A"], pops["C"], pops["M"], pops["O"]], indexes=[(0,1), (0,2), (0, 3), (1, 2), (1, 3), (2,3)])
#print(f2)
f2DF = pd.DataFrame(f2, columns = ["f2"])
f2DF.loc[:, "Group1"] = ["A", "A", "A", "C", "C", "M"]
f2DF.loc[:, "Group2"] = ["C", "M", "O", "M", "O", "O"]
f2DF.to_csv("f2.csv", index=None)

# f4
f4 = drone_ts.f4([pops["A"], pops["C"], pops["M"], pops["O"]])
print("This is f4 statistics: ", str(f4))


# In[ ]:


# GNN
gnn = pd.DataFrame(drone_ts.genealogical_nearest_neighbours(pops["H"], sample_sets=[pops["A"], pops["C"], pops["M"], pops["O"]]), columns = ["A", "C","M", "O"])

# Add sample ids and names to the data frame
#gnn.loc[:, "Id"] = pops["H"] #The two nodes per individual will be the same since they are homozygous diploids
namesH = [name for (name, ID) in zip(sample_names, sample_ids) if ID in pops['H']]
gnn.loc[:, "Names"] = namesH
gnn = gnn.drop_duplicates()

gnn.to_csv("GNN_drones.csv", index=None)


# In[ ]:


# Genetic relatedness
grel = drone_ts.genetic_relatedness([pops["A"], pops["C"], pops["M"], pops["O"]], indexes=[(0,1), (0,2), (0, 3), (1, 2), (1, 3), (2,3)])
#print(grel)
grelDF = pd.DataFrame(grel, columns = ["GeneticRel"])
grelDF.loc[:, "Group1"] = ["A", "A", "A", "C", "C", "M"]
grelDF.loc[:, "Group2"] = ["C", "M", "O", "M", "O", "O"]
grelDF.to_csv("GeneticRelatedness.csv", index=None)


# In[ ]:


# Mean descendants of ALL nodes!!! Also intermediate
meanDesc = pd.DataFrame(drone_ts.mean_descendants([pops["A"], pops["C"], pops["M"], pops["O"]]), columns = ["A", "C", "M", "O"])
meanDesc.loc[:, "NodeID"] = [n.id for n in drone_ts.nodes()]
meanDesc.loc[:, "NodeName"] = None

# Add names for sample nodes
for (ID, name) in zip(sample_nodes, sample_names):
    meanDesc.at[meanDesc.NodeID == ID.id, "NodeName"] = name
    

meanDesc.to_csv("MeanDescendants.csv", index=None)
            
            
        


# In[ ]:


# Write the tables with mutations, nodes, sites, populations and individuals
drone_ts.dump_text(mutations=open("Mutations.txt", "w"), nodes = open("Nodes.txt", "w"), sites = open("Sites.txt", "w"), edges = open("Edges.txt", "w"), individuals = open("Individuals.txt", "w"), populations = open("Populations.txt", "w", encoding="utf-8"))

