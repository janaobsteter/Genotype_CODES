SELEKCIJA_GUI.py

1) SELECTION PARAMETERS
No Newborns - število novorojenih ženskih in moških telet v eni generaciji

ŽENSKA LINIJA
% of NB > female calves: kakšen % novorojenih Ž živali preživi do prvega leta (so telice)
% females calves > cows: kakšen % telic osemenimo
% dams (cows): kakšen % krav postane bikovskih mater

MOŠKA LINIJA
% načrtnih parjenj (od NB): kakšen % novorojenih moških je potomcev načrtnih parjenj
% classical test / genomic testing: ali imajo bikci progeni test (potomci NP > vhlevljeni > mladi > pozitivno testirani) 
ali so genomsko testirani (vsi potomci NP gEBV > pozitivno testirani)
	ČE PROGENI TEST
	% vhlevljeni: kakšen % potomcev NP je vhlevljenih
	% mladi: kakšen % vhlevljnih bikov postane mladih (XX doz)
% pozTest: kakešn % mladih / genomsko testiranih bikov prestane test in se uporablja v AI
% pripust: kakšen % bikov gre za pripust (tu mora biti mladi + pripust = 100%)
% of NB > male calves: kakšen % novorojenih moških živali preživi do prvega leta
% of male calves > bulls: kakšen % novorojenih moških živali preživi do drugega leta

NATANČNOSTI
Accuracy of EBV / gEBV: natančnost napovedi PV

2) TIME PARAMETERS
Lactations per cow: povprečno koliko laktacij ima krava
Dams chosen after lactaion no: po kateri lactaciji odbiramo bikovske matere
Lactations per dam: povprečno koliko laktacij so krave bikovske matere
years in test (bulls): koliko časa so biki v testu (progeni test) - ne štejemo leta, ko so mladi bikI!
natural service in use: koliko let so v uporabi biki v pripustu
testes bulls years in use: koliko let so v uporabi pozitivno testirani biki

3) DOSAGE PARAMETERS
Dosage per young bulls: koliko doz po mladem biku za AI (progeni test)
Dosage per natural service bull: koliko doz po biku v pripustu
Dosage per tested: koliko doz po pozitivno testiranem biku za AI (progeni test ali genomsko testiran)

4) ALPHASIM
PErform burn in: ali že imaš burn in populacijo ali jo naredi
	NoSires, NoDams: koliko je mater / očetov v eni generaciji
Plot genetic gain: da/ne proizvedi graf povprečnih TBV skozi generacije
NoBurnIn: koliko je burn in generacij (vnesti tudi, če je burn in že narejen)
No Selected Generation: koliko generacij želiš izvajati selekcijo
Perform selection From/To Generation: koliko in katere cikle selekcije hočeš pognati
AlphaSim: najdi direktorij, kjer se nahaja AlphaSim
