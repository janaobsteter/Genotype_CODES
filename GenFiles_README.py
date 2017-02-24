GenFiles.py module documentation

Dictionaries: 
-dictionary of numSNP: chip
-TraitSNPs = {Trait: {Subtrait : SNPs}}
-SNP800Sifrant_Dict = {SNP : {sifra1, sifra2}}


classes:
genZipPackage: 
-__init__ = self.name, self.genodate, self.infiles (lists files in the zip)
- unzip() - extracts all from the zip file
- extractFinalReport - extracts FinalReport, if only one the package, no subdirectories will be created for the regarding genotyping package. If multiple genotype packages within, subdirectory for each will be created
- extractSNPMap: same as for the FinalReport
- checkSubDir: check, whether the zip file includes multiple subdirectories, which will be created at the extraction of the files

pedFile:
-__init__ = name, sernum, genodate, pedContent (reads the files), samples (sample IDs), chip, mapContent, snps (lists the SNPs from the corresponding map file)
-extractSNP(SNP) = extract one SNP
-extractSNPList (SNPlist = list) = extract a list of chosen SNPs, enter them as a list
-extract TratSNPs (Trait) = specify the Trait you want to extract SNPs for (source of SNPs is the dictionary)
-extractParentalSNPs(number): extracts the chosen palet of SNPs for parental verification

mapFiles:
-__init__ = name, sernum, mapContent, snps, chip
-chrSNP(chromosome): count and lists (optional) the SNPs on a specific chromosome
-posSNP(chromosome, pos): lists SNPs on the specific position
