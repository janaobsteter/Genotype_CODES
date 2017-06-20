# -*- coding: utf-8 -*-
blupFiles = blupf90(AlphaSimDir, way='burnin_milk', sel='gen')
blupFiles.makeDat_sex(2)
shutil.copy(blupFiles.blupgenParamFile, blupFiles.AlphaSimDir)  # skopiraj template blupparam file

# uredi blupparam file
# get variance components from AlphaSim Output Files
OutputFiles = AlphaSim_OutputFile(AlphaSimDir)
genvar = OutputFiles.getAddVar()  # dobi additivno varianco
resvar = OutputFiles.getResVar()  # dobi varianco za ostanek

blupFiles.prepareParamFiles(genvar, resvar)  # set levels of random aniaml effect, add var and res var
# the paramfile is now set
if sel == 'class':
    blupFiles.makePed_class()  # make ped file for blup, code (1, 2, 3 - both parentr knows/unknown/group)
    os.system('./blupf90 blupf90_Selection')
if sel == 'gen':
    blupFiles.makePed_gen()  # make ped file for blup, no Code!
    GenFiles = snpFiles(AlphaSimDir)
    GenFiles.createBlupf90SNPFile()
    os.system('./renumf90 < renumParam')  # run blupf90
    # os.system('./blupf90 blupf90_Selection')
    resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    os.system('./preGSf90 renf90.par')
    os.system('./blupf90 renf90.par')
    # os.system('./postGSf90 renf90.par')
    # os.system('./ renf90.par')
    
    
    
blupSol = pd.read_csv(AlphaSimDir + '/solutions', skiprows=1, header=None, sep='\s+', names=['Trait', 'Effect', 'Level', 'Solution'])
blupSol = pd.read_csv(AlphaSimDir + '/renumbered_Solutions', header=None, sep='\s+', names=['RenumID', 'OrigID', 'Solution'])



blupSol = pd.read_csv(AlphaSimDir + '/renumbered_Solutions', skiprows=1, header=None, 
                        sep='\s+', names=['renID', 'ID', 'Solution'])
AlphaSelPed = AlphaPed.loc[:, ['Generation', 'Indiv', 'Father', 'Mother', 'gvNormUnres1']]
#blupSolRandom = blupSol.loc[blupSol.Effect == 1] Če imaš še fixed effect
AlphaSelPed.loc[:, 'EBV'] = blupSol.Solution
AlphaSelPed.to_csv(AlphaSimDir + 'GenPed_EBV.txt', index=None)



if os.path.isfile(AlphaSimDir + 'GenoFile.txt'):
    os.system('''sed -n "$(sed 's/$/p/' IndForGeno_new.txt)" ''' + chip + ' > ChosenInd.txt') #only individuals chosen for genotypisation - ONLY NEW - LAST GEN!
    os.system("sed 's/^ *//' ChosenInd.txt > ChipFile.txt") #Remove blank spaces at the beginning
    os.system("cut -f1 -d ' ' ChipFile.txt > Individuals.txt")  # obtain IDs
    os.system('''awk '{$1=""; print $0}' ChipFile.txt | sed 's/ //g' > Snps.txt''')  # obtain SNP genotypes
    os.system(
        r'''paste Individuals.txt Snps.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' > GenoFile_new.txt''')  # obtain SNP genotypes of the last generation
    os.system("cat GenoFile.txt GenoFile_new.txt > GenoFileTmp && mv GenoFileTmp GenoFile.txt")
else:
    os.system('''sed -n "$(sed 's/$/p/' IndForGeno.txt)" ''' + chip + ' > ChosenInd.txt') #only individuals chosen for genotypisation - ALL
    os.system("sed 's/^ *//' ChosenInd.txt > ChipFile.txt") #Remove blank spaces at the beginning
    os.system("cut -f1 -d ' ' ChipFile.txt > Individuals.txt") #obtain IDs
    os.system('''awk '{$1=""; print $0}' ChipFile.txt | sed 's/ //g' > Snps.txt''') #obtain SNP genotypes
    os.system(
        r'''paste Individuals.txt Snps.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' > GenoFile.txt''')  # obtain SNP genotypes of the last generation
pd.read_csv(AlphaSimDir + '/SimulatedData/Chip1SnpInformation.txt', sep='\s+')[[0, 1, 2]].to_csv(AlphaSimDir + 'SnpMap.txt', index=None, sep=" ", header=None)
