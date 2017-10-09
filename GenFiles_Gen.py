# -*- coding: utf-8 -*-
import os
from collections import defaultdict
import csv
import pandas as pd
import zipfile
import shutil


chips = {19720: "GGPv02", 
26145: "GGPv03", 
26151: "GGPv03", 
30105: "GGPv04",
30106: "GGPv04", 
76883: "HD" , 
138892: "HDv02", 
139376: "HDv02", 
54001:"50Kv01" , 
54609: "50Kv02",
51274: "IDBv03"
         }

TraitSNPs = {
'Caseins':['BCNAB', 'BCNAB_2', 'BCNAB_3','GNSC319','GNSC319_3','GNSC319_B1','GNSC355','GNSC355_3','GNSC355_B1','KappaCasein12951_1','KappaCasein12951_2','KappaCasein12951_3'],
'BetaCasein' : ['BCNAB', 'BCNAB_2', 'BCNAB_3'], 
'KappaCasein': ['GNSC319','GNSC319_3','GNSC319_B1','GNSC355','GNSC355_3','GNSC355_B1','KappaCasein12951_1','KappaCasein12951_2','KappaCasein12951_3'],
'Citrulinemia':[]}

Parental800=open("/home/jana/Genotipi/ParentalVerification_SNPSNP/Names_800SNPs.txt").read().strip("\n").split("\n")



"""
To je samo za ustvarit Sifrant SNPov,potem, ko ga imaš, ga skopiraj v datoteko
SNPSifrant="/home/jana/Genotipi/ParentalVerification_SNPSNP/Sifrant_SNP.csv"

SNPSifrant_Dict = defaultdict(set)
with open(SNPSifrant, 'rb') as SNP_Sifrant:
    SNP_Sifrant.readline() #skip the first line
    reader = csv.reader(SNP_Sifrant, delimiter=',')
    for line in reader:
        SNPSifrant_Dict[line[0]].add(line[1])
        
"""
SNP800Sifrant_Dict = {
'ARS-BFGL-BAC-11276': {'1', '2'},
'ARS-BFGL-BAC-12893': {'3', '4'},
'ARS-BFGL-BAC-14220': {'5', '6'},
'ARS-BFGL-BAC-14359': {'7', '8'},
'ARS-BFGL-BAC-14813': {'10', '9'},
'ARS-BFGL-BAC-15437': {'11', '12'},
'ARS-BFGL-BAC-16276': {'13', '14'},
'ARS-BFGL-BAC-19454': {'15', '16'},
'ARS-BFGL-BAC-23417': {'17', '18'},
'ARS-BFGL-BAC-2599': {'19', '20'},
'ARS-BFGL-BAC-27364': {'21', '22'},
'ARS-BFGL-BAC-29664': {'23', '24'},
'ARS-BFGL-BAC-30348': {'25', '26'},
'ARS-BFGL-BAC-30716': {'27', '28'},
'ARS-BFGL-BAC-32404': {'29', '30'},
'ARS-BFGL-BAC-35552': {'31', '32'},
'ARS-BFGL-BAC-36979': {'33', '34'},
'ARS-BFGL-BAC-40619': {'35', '36'},
'ARS-BFGL-BAC-6489': {'37', '38'},
'ARS-BFGL-NGS-100131': {'39', '40'},
'ARS-BFGL-NGS-10035': {'41', '42'},
'ARS-BFGL-NGS-10128': {'43', '44'},
'ARS-BFGL-NGS-101456': {'45', '46'},
'ARS-BFGL-NGS-101468': {'47', '48'},
'ARS-BFGL-NGS-101742': {'49', '50'},
'ARS-BFGL-NGS-101968': {'51', '52'},
'ARS-BFGL-NGS-102169': {'53', '54'},
'ARS-BFGL-NGS-10254': {'55', '56'},
'ARS-BFGL-NGS-102603': {'57', '58'},
'ARS-BFGL-NGS-102610': {'59', '60'},
'ARS-BFGL-NGS-102877': {'61', '62'},
'ARS-BFGL-NGS-103280': {'63', '64'},
'ARS-BFGL-NGS-103379': {'65', '66'},
'ARS-BFGL-NGS-103387': {'67', '68'},
'ARS-BFGL-NGS-103449': {'69', '70'},
'ARS-BFGL-NGS-103469': {'71', '72'},
'ARS-BFGL-NGS-103524': {'73', '74'},
'ARS-BFGL-NGS-103691': {'75', '76'},
'ARS-BFGL-NGS-103861': {'77', '78'},
'ARS-BFGL-NGS-10436': {'79', '80'},
'ARS-BFGL-NGS-104603': {'81', '82'},
'ARS-BFGL-NGS-104897': {'83', '84'},
'ARS-BFGL-NGS-10560': {'85', '86'},
'ARS-BFGL-NGS-105811': {'87', '88'},
'ARS-BFGL-NGS-106015': {'89', '90'},
'ARS-BFGL-NGS-106289': {'91', '92'},
'ARS-BFGL-NGS-106379': {'93', '94'},
'ARS-BFGL-NGS-107085': {'95', '96'},
'ARS-BFGL-NGS-107176': {'97', '98'},
'ARS-BFGL-NGS-107540': {'100', '99'},
'ARS-BFGL-NGS-107714': {'101', '102'},
'ARS-BFGL-NGS-107767': {'103', '104'},
'ARS-BFGL-NGS-108099': {'105', '106'},
'ARS-BFGL-NGS-108308': {'107', '108'},
'ARS-BFGL-NGS-108549': {'109', '110'},
'ARS-BFGL-NGS-108579': {'111', '112'},
'ARS-BFGL-NGS-109532': {'113', '114'},
'ARS-BFGL-NGS-1097': {'115', '116'},
'ARS-BFGL-NGS-109722': {'117', '118'},
'ARS-BFGL-NGS-109750': {'119', '120'},
'ARS-BFGL-NGS-109752': {'121', '122'},
'ARS-BFGL-NGS-109907': {'123', '124'},
'ARS-BFGL-NGS-110190': {'125', '126'},
'ARS-BFGL-NGS-110280': {'127', '128'},
'ARS-BFGL-NGS-110420': {'129', '130'},
'ARS-BFGL-NGS-110438': {'131', '132'},
'ARS-BFGL-NGS-110665': {'133', '134'},
'ARS-BFGL-NGS-110811': {'135', '136'},
'ARS-BFGL-NGS-110930': {'137', '138'},
'ARS-BFGL-NGS-110965': {'139', '140'},
'ARS-BFGL-NGS-111049': {'141', '142'},
'ARS-BFGL-NGS-111053': {'143', '144'},
'ARS-BFGL-NGS-111076': {'145', '146'},
'ARS-BFGL-NGS-111114': {'147', '148'},
'ARS-BFGL-NGS-111118': {'149', '150'},
'ARS-BFGL-NGS-111681': {'151', '152'},
'ARS-BFGL-NGS-111691': {'153', '154'},
'ARS-BFGL-NGS-111777': {'155', '156'},
'ARS-BFGL-NGS-111788': {'157', '158'},
'ARS-BFGL-NGS-111859': {'159', '160'},
'ARS-BFGL-NGS-112094': {'161', '162'},
'ARS-BFGL-NGS-112130': {'163', '164'},
'ARS-BFGL-NGS-112133': {'165', '166'},
'ARS-BFGL-NGS-112325': {'167', '168'},
'ARS-BFGL-NGS-112372': {'169', '170'},
'ARS-BFGL-NGS-112392': {'171', '172'},
'ARS-BFGL-NGS-112592': {'173', '174'},
'ARS-BFGL-NGS-112625': {'175', '176'},
'ARS-BFGL-NGS-113524': {'177', '178'},
'ARS-BFGL-NGS-113570': {'179', '180'},
'ARS-BFGL-NGS-113716': {'181', '182'},
'ARS-BFGL-NGS-11383': {'183', '184'},
'ARS-BFGL-NGS-113910': {'185', '186'},
'ARS-BFGL-NGS-113972': {'187', '188'},
'ARS-BFGL-NGS-114006': {'189', '190'},
'ARS-BFGL-NGS-114028': {'191', '192'},
'ARS-BFGL-NGS-114067': {'193', '194'},
'ARS-BFGL-NGS-114087': {'195', '196'},
'ARS-BFGL-NGS-114245': {'197', '198'},
'ARS-BFGL-NGS-114538': {'199', '200'},
'ARS-BFGL-NGS-114766': {'201', '202'},
'ARS-BFGL-NGS-114947': {'203', '204'},
'ARS-BFGL-NGS-115237': {'205', '206'},
'ARS-BFGL-NGS-115346': {'207', '208'},
'ARS-BFGL-NGS-115377': {'209', '210'},
'ARS-BFGL-NGS-115514': {'211', '212'},
'ARS-BFGL-NGS-115598': {'213', '214'},
'ARS-BFGL-NGS-115613': {'215', '216'},
'ARS-BFGL-NGS-11585': {'217', '218'},
'ARS-BFGL-NGS-116095': {'219', '220'},
'ARS-BFGL-NGS-116165': {'221', '222'},
'ARS-BFGL-NGS-116386': {'223', '224'},
'ARS-BFGL-NGS-116541': {'225', '226'},
'ARS-BFGL-NGS-116543': {'227', '228'},
'ARS-BFGL-NGS-116589': {'229', '230'},
'ARS-BFGL-NGS-116761': {'231', '232'},
'ARS-BFGL-NGS-117122': {'233', '234'},
'ARS-BFGL-NGS-117241': {'235', '236'},
'ARS-BFGL-NGS-117319': {'237', '238'},
'ARS-BFGL-NGS-117322': {'239', '240'},
'ARS-BFGL-NGS-117421': {'241', '242'},
'ARS-BFGL-NGS-117553': {'243', '244'},
'ARS-BFGL-NGS-117619': {'245', '246'},
'ARS-BFGL-NGS-117672': {'247', '248'},
'ARS-BFGL-NGS-117870': {'249', '250'},
'ARS-BFGL-NGS-118073': {'251', '252'},
'ARS-BFGL-NGS-118102': {'253', '254'},
'ARS-BFGL-NGS-118138': {'255', '256'},
'ARS-BFGL-NGS-118319': {'257', '258'},
'ARS-BFGL-NGS-118340': {'259', '260'},
'ARS-BFGL-NGS-118398': {'261', '262'},
'ARS-BFGL-NGS-118453': {'263', '264'},
'ARS-BFGL-NGS-118597': {'265', '266'},
'ARS-BFGL-NGS-118621': {'267', '268'},
'ARS-BFGL-NGS-118876': {'269', '270'},
'ARS-BFGL-NGS-11889': {'271', '272'},
'ARS-BFGL-NGS-119431': {'273', '274'},
'ARS-BFGL-NGS-119662': {'275', '276'},
'ARS-BFGL-NGS-119694': {'277', '278'},
'ARS-BFGL-NGS-119782': {'279', '280'},
'ARS-BFGL-NGS-12481': {'281', '282'},
'ARS-BFGL-NGS-12514': {'283', '284'},
'ARS-BFGL-NGS-12595': {'285', '286'},
'ARS-BFGL-NGS-12759': {'287', '288'},
'ARS-BFGL-NGS-12878': {'289', '290'},
'ARS-BFGL-NGS-13360': {'291', '292'},
'ARS-BFGL-NGS-14035': {'293', '294'},
'ARS-BFGL-NGS-14039': {'295', '296'},
'ARS-BFGL-NGS-14557': {'297', '298'},
'ARS-BFGL-NGS-14632': {'299', '300'},
'ARS-BFGL-NGS-14740': {'301', '302'},
'ARS-BFGL-NGS-1481': {'303', '304'},
'ARS-BFGL-NGS-14996': {'305', '306'},
'ARS-BFGL-NGS-15250': {'307', '308'},
'ARS-BFGL-NGS-15506': {'309', '310'},
'ARS-BFGL-NGS-15508': {'311', '312'},
'ARS-BFGL-NGS-15731': {'313', '314'},
'ARS-BFGL-NGS-16249': {'315', '316'},
'ARS-BFGL-NGS-16336': {'317', '318'},
'ARS-BFGL-NGS-16504': {'319', '320'},
'ARS-BFGL-NGS-16920': {'321', '322'},
'ARS-BFGL-NGS-16925': {'323', '324'},
'ARS-BFGL-NGS-17077': {'325', '326'},
'ARS-BFGL-NGS-17125': {'327', '328'},
'ARS-BFGL-NGS-17182': {'329', '330'},
'ARS-BFGL-NGS-17527': {'331', '332'},
'ARS-BFGL-NGS-18122': {'333', '334'},
'ARS-BFGL-NGS-18653': {'335', '336'},
'ARS-BFGL-NGS-19227': {'337', '338'},
'ARS-BFGL-NGS-19431': {'339', '340'},
'ARS-BFGL-NGS-19741': {'341', '342'},
'ARS-BFGL-NGS-20281': {'343', '344'},
'ARS-BFGL-NGS-20403': {'345', '346'},
'ARS-BFGL-NGS-21227': {'347', '348'},
'ARS-BFGL-NGS-21400': {'349', '350'},
'ARS-BFGL-NGS-21717': {'351', '352'},
'ARS-BFGL-NGS-22052': {'353', '354'},
'ARS-BFGL-NGS-22113': {'355', '356'},
'ARS-BFGL-NGS-2215': {'357', '358'},
'ARS-BFGL-NGS-22685': {'359', '360'},
'ARS-BFGL-NGS-23432': {'361', '362'},
'ARS-BFGL-NGS-23813': {'363', '364'},
'ARS-BFGL-NGS-23961': {'365', '366'},
'ARS-BFGL-NGS-24214': {'367', '368'},
'ARS-BFGL-NGS-24419': {'369', '370'},
'ARS-BFGL-NGS-24521': {'371', '372'},
'ARS-BFGL-NGS-24591': {'373', '374'},
'ARS-BFGL-NGS-25107': {'375', '376'},
'ARS-BFGL-NGS-25646': {'377', '378'},
'ARS-BFGL-NGS-25959': {'379', '380'},
'ARS-BFGL-NGS-2644': {'381', '382'},
'ARS-BFGL-NGS-26517': {'383', '384'},
'ARS-BFGL-NGS-27095': {'385', '386'},
'ARS-BFGL-NGS-27577': {'387', '388'},
'ARS-BFGL-NGS-27731': {'389', '390'},
'ARS-BFGL-NGS-27865': {'391', '392'},
'ARS-BFGL-NGS-27948': {'393', '394'},
'ARS-BFGL-NGS-28167': {'395', '396'},
'ARS-BFGL-NGS-28202': {'397', '398'},
'ARS-BFGL-NGS-28504': {'399', '400'},
'ARS-BFGL-NGS-28681': {'401', '402'},
'ARS-BFGL-NGS-29144': {'403', '404'},
'ARS-BFGL-NGS-30199': {'405', '406'},
'ARS-BFGL-NGS-30207': {'407', '408'},
'ARS-BFGL-NGS-30559': {'409', '410'},
'ARS-BFGL-NGS-30677': {'411', '412'},
'ARS-BFGL-NGS-30953': {'413', '414'},
'ARS-BFGL-NGS-31640': {'415', '416'},
'ARS-BFGL-NGS-31807': {'417', '418'},
'ARS-BFGL-NGS-32083': {'419', '420'},
'ARS-BFGL-NGS-32490': {'421', '422'},
'ARS-BFGL-NGS-32605': {'423', '424'},
'ARS-BFGL-NGS-32769': {'425', '426'},
'ARS-BFGL-NGS-32846': {'427', '428'},
'ARS-BFGL-NGS-33000': {'429', '430'},
'ARS-BFGL-NGS-33249': {'431', '432'},
'ARS-BFGL-NGS-33709': {'433', '434'},
'ARS-BFGL-NGS-33745': {'435', '436'},
'ARS-BFGL-NGS-34169': {'437', '438'},
'ARS-BFGL-NGS-34282': {'439', '440'},
'ARS-BFGL-NGS-34393': {'441', '442'},
'ARS-BFGL-NGS-34516': {'443', '444'},
'ARS-BFGL-NGS-34582': {'445', '446'},
'ARS-BFGL-NGS-35419': {'447', '448'},
'ARS-BFGL-NGS-3584': {'449', '450'},
'ARS-BFGL-NGS-35963': {'451', '452'},
'ARS-BFGL-NGS-36513': {'453', '454'},
'ARS-BFGL-NGS-36545': {'455', '456'},
'ARS-BFGL-NGS-36553': {'457', '458'},
'ARS-BFGL-NGS-36793': {'459', '460'},
'ARS-BFGL-NGS-36852': {'461', '462'},
'ARS-BFGL-NGS-36880': {'463', '464'},
'ARS-BFGL-NGS-37484': {'465', '466'},
'ARS-BFGL-NGS-37680': {'467', '468'},
'ARS-BFGL-NGS-38423': {'469', '470'},
'ARS-BFGL-NGS-38620': {'471', '472'},
'ARS-BFGL-NGS-3965': {'473', '474'},
'ARS-BFGL-NGS-39913': {'475', '476'},
'ARS-BFGL-NGS-39978': {'477', '478'},
'ARS-BFGL-NGS-40189': {'479', '480'},
'ARS-BFGL-NGS-40229': {'481', '482'},
'ARS-BFGL-NGS-40249': {'483', '484'},
'ARS-BFGL-NGS-4047': {'485', '486'},
'ARS-BFGL-NGS-40856': {'487', '488'},
'ARS-BFGL-NGS-41117': {'489', '490'},
'ARS-BFGL-NGS-41462': {'491', '492'},
'ARS-BFGL-NGS-42283': {'493', '494'},
'ARS-BFGL-NGS-42338': {'495', '496'},
'ARS-BFGL-NGS-42505': {'497', '498'},
'ARS-BFGL-NGS-42852': {'499', '500'},
'ARS-BFGL-NGS-4292': {'501', '502'},
'ARS-BFGL-NGS-42950': {'503', '504'},
'ARS-BFGL-NGS-44153': {'505', '506'},
'ARS-BFGL-NGS-44355': {'507', '508'},
'ARS-BFGL-NGS-44611': {'509', '510'},
'ARS-BFGL-NGS-44686': {'511', '512'},
'ARS-BFGL-NGS-45453': {'513', '514'},
'ARS-BFGL-NGS-45547': {'515', '516'},
'ARS-BFGL-NGS-45686': {'517', '518'},
'ARS-BFGL-NGS-50172': {'519', '520'},
'ARS-BFGL-NGS-5022': {'521', '522'},
'ARS-BFGL-NGS-5133': {'523', '524'},
'ARS-BFGL-NGS-53433': {'525', '526'},
'ARS-BFGL-NGS-53461': {'527', '528'},
'ARS-BFGL-NGS-54753': {'529', '530'},
'ARS-BFGL-NGS-55042': {'531', '532'},
'ARS-BFGL-NGS-55943': {'533', '534'},
'ARS-BFGL-NGS-56579': {'535', '536'},
'ARS-BFGL-NGS-5662': {'537', '538'},
'ARS-BFGL-NGS-5665': {'539', '540'},
'ARS-BFGL-NGS-57149': {'541', '542'},
'ARS-BFGL-NGS-57417': {'543', '544'},
'ARS-BFGL-NGS-57549': {'545', '546'},
'ARS-BFGL-NGS-57711': {'547', '548'},
'ARS-BFGL-NGS-57870': {'549', '550'},
'ARS-BFGL-NGS-58066': {'551', '552'},
'ARS-BFGL-NGS-58613': {'553', '554'},
'ARS-BFGL-NGS-60640': {'555', '556'},
'ARS-BFGL-NGS-6240': {'557', '558'},
'ARS-BFGL-NGS-6269': {'559', '560'},
'ARS-BFGL-NGS-63132': {'561', '562'},
'ARS-BFGL-NGS-64451': {'563', '564'},
'ARS-BFGL-NGS-64841': {'565', '566'},
'ARS-BFGL-NGS-66512': {'567', '568'},
'ARS-BFGL-NGS-67112': {'569', '570'},
'ARS-BFGL-NGS-67146': {'571', '572'},
'ARS-BFGL-NGS-67161': {'573', '574'},
'ARS-BFGL-NGS-67260': {'575', '576'},
'ARS-BFGL-NGS-67754': {'577', '578'},
'ARS-BFGL-NGS-67968': {'579', '580'},
'ARS-BFGL-NGS-68030': {'581', '582'},
'ARS-BFGL-NGS-68850': {'583', '584'},
'ARS-BFGL-NGS-70946': {'585', '586'},
'ARS-BFGL-NGS-71395': {'587', '588'},
'ARS-BFGL-NGS-72192': {'589', '590'},
'ARS-BFGL-NGS-72413': {'591', '592'},
'ARS-BFGL-NGS-72471': {'593', '594'},
'ARS-BFGL-NGS-73980': {'595', '596'},
'ARS-BFGL-NGS-74123': {'597', '598'},
'ARS-BFGL-NGS-7433': {'599', '600'},
'ARS-BFGL-NGS-74380': {'601', '602'},
'ARS-BFGL-NGS-74523': {'603', '604'},
'ARS-BFGL-NGS-74971': {'605', '606'},
'ARS-BFGL-NGS-75066': {'607', '608'},
'ARS-BFGL-NGS-75391': {'609', '610'},
'ARS-BFGL-NGS-75852': {'611', '612'},
'ARS-BFGL-NGS-76111': {'613', '614'},
'ARS-BFGL-NGS-76330': {'615', '616'},
'ARS-BFGL-NGS-80026': {'617', '618'},
'ARS-BFGL-NGS-80841': {'619', '620'},
'ARS-BFGL-NGS-81155': {'621', '622'},
'ARS-BFGL-NGS-81428': {'623', '624'},
'ARS-BFGL-NGS-8229': {'625', '626'},
'ARS-BFGL-NGS-82717': {'627', '628'},
'ARS-BFGL-NGS-83035': {'629', '630'},
'ARS-BFGL-NGS-83273': {'631', '632'},
'ARS-BFGL-NGS-84065': {'633', '634'},
'ARS-BFGL-NGS-84327': {'635', '636'},
'ARS-BFGL-NGS-8471': {'637', '638'},
'ARS-BFGL-NGS-86662': {'639', '640'},
'ARS-BFGL-NGS-87679': {'641', '642'},
'ARS-BFGL-NGS-89098': {'643', '644'},
'ARS-BFGL-NGS-89746': {'645', '646'},
'ARS-BFGL-NGS-8982': {'647', '648'},
'ARS-BFGL-NGS-92824': {'649', '650'},
'ARS-BFGL-NGS-92946': {'651', '652'},
'ARS-BFGL-NGS-93119': {'653', '654'},
'ARS-BFGL-NGS-94543': {'655', '656'},
'ARS-BFGL-NGS-94895': {'657', '658'},
'ARS-BFGL-NGS-95360': {'659', '660'},
'ARS-BFGL-NGS-96125': {'661', '662'},
'ARS-BFGL-NGS-98526': {'663', '664'},
'ARS-BFGL-NGS-99210': {'665', '666'},
'ARS-USMARC-569': {'667', '668'},
'ARS-USMARC-Parent-AY761135-rs29003723': {'669', '670'},
'ARS-USMARC-Parent-AY776154-no-rs': {'671', '672'},
'ARS-USMARC-Parent-AY841151-rs29003466': {'673', '674'},
'ARS-USMARC-Parent-AY842472-rs29001941': {'675', '676'},
'ARS-USMARC-Parent-AY842473-rs29001956': {'677', '678'},
'ARS-USMARC-Parent-AY842474-rs29003226': {'679', '680'},
'ARS-USMARC-Parent-AY842475-rs29002127': {'681', '682'},
'ARS-USMARC-Parent-AY844963-rs17871338': {'683', '684'},
'ARS-USMARC-Parent-AY849381-rs29003287': {'685', '686'},
'ARS-USMARC-Parent-AY850194-no-rs': {'687', '688'},
'ARS-USMARC-Parent-AY851162-no-rs': {'689', '690'},
'ARS-USMARC-Parent-AY851163-rs17871661': {'691', '692'},
'ARS-USMARC-Parent-AY853302-no-rs': {'693', '694'},
'ARS-USMARC-Parent-AY853303-no-rs': {'695', '696'},
'ARS-USMARC-Parent-AY856094-rs17871190': {'697', '698'},
'ARS-USMARC-Parent-AY857620-rs17871214': {'699', '700'},
'ARS-USMARC-Parent-AY858890-rs29002256': {'701', '702'},
'ARS-USMARC-Parent-AY860426-no-rs': {'703', '704'},
'ARS-USMARC-Parent-AY863214-rs17871744': {'705', '706'},
'ARS-USMARC-Parent-AY914316-rs17871403': {'707', '708'},
'ARS-USMARC-Parent-AY916666-no-rs': {'709', '710'},
'ARS-USMARC-Parent-AY919868-rs29002211': {'711', '712'},
'ARS-USMARC-Parent-AY929334-no-rs': {'713', '714'},
'ARS-USMARC-Parent-AY937242-rs17872223': {'715', '716'},
'ARS-USMARC-Parent-AY939849-rs17870274': {'717', '718'},
'ARS-USMARC-Parent-AY941204-rs17872131': {'719', '720'},
'ARS-USMARC-Parent-AY943841-rs17871566': {'721', '722'},
'ARS-USMARC-Parent-DQ381152-rs29002408': {'723', '724'},
'ARS-USMARC-Parent-DQ381153-rs29012842': {'725', '726'},
'ARS-USMARC-Parent-DQ404149-no-rs': {'727', '728'},
'ARS-USMARC-Parent-DQ404150-rs29012530': {'729', '730'},
'ARS-USMARC-Parent-DQ404151-rs29019282': {'731', '732'},
'ARS-USMARC-Parent-DQ404152-rs29022245': {'733', '734'},
'ARS-USMARC-Parent-DQ404153-no-rs': {'735', '736'},
'ARS-USMARC-Parent-DQ435443-rs29010802': {'737', '738'},
'ARS-USMARC-Parent-DQ451555-rs29010795': {'739', '740'},
'ARS-USMARC-Parent-DQ468384-rs29003967': {'741', '742'},
'ARS-USMARC-Parent-DQ470475-no-rs': {'743', '744'},
'ARS-USMARC-Parent-DQ489377-rs29026932': {'745', '746'},
'ARS-USMARC-Parent-DQ500958-no-rs': {'747', '748'},
'ARS-USMARC-Parent-DQ647186-rs29014143': {'749', '750'},
'ARS-USMARC-Parent-DQ647187-rs29010510': {'751', '752'},
'ARS-USMARC-Parent-DQ647188-rs29011099': {'753', '754'},
'ARS-USMARC-Parent-DQ647189-rs29012226': {'755', '756'},
'ARS-USMARC-Parent-DQ647190-rs29013632': {'757', '758'},
'ARS-USMARC-Parent-DQ650635-rs29012174': {'759', '760'},
'ARS-USMARC-Parent-DQ650636-rs29024525': {'761', '762'},
'ARS-USMARC-Parent-DQ674265-rs29011266': {'763', '764'},
'ARS-USMARC-Parent-DQ786757-rs29019900': {'765', '766'},
'ARS-USMARC-Parent-DQ786758-rs29024430': {'767', '768'},
'ARS-USMARC-Parent-DQ786759-rs29026696': {'769', '770'},
'ARS-USMARC-Parent-DQ786761-rs29012840': {'771', '772'},
'ARS-USMARC-Parent-DQ786762-rs29010772': {'773', '774'},
'ARS-USMARC-Parent-DQ786763-rs29020472': {'775', '776'},
'ARS-USMARC-Parent-DQ786765-rs29009858': {'777', '778'},
'ARS-USMARC-Parent-DQ786766-rs29012070': {'779', '780'},
'ARS-USMARC-Parent-DQ789028-rs29017713': {'781', '782'},
'ARS-USMARC-Parent-DQ837643-rs29018818': {'783', '784'},
'ARS-USMARC-Parent-DQ837644-rs29010468': {'785', '786'},
'ARS-USMARC-Parent-DQ839235-rs29012691': {'787', '788'},
'ARS-USMARC-Parent-DQ846688-rs29023691': {'789', '790'},
'ARS-USMARC-Parent-DQ846690-no-rs': {'791', '792'},
'ARS-USMARC-Parent-DQ846691-rs29019814': {'793', '794'},
'ARS-USMARC-Parent-DQ846692-rs29010281': {'795', '796'},
'ARS-USMARC-Parent-DQ846693-rs29017621': {'797', '798'},
'ARS-USMARC-Parent-DQ866817-no-rs': {'799', '800'},
'ARS-USMARC-Parent-DQ866818-rs29011701': {'801', '802'},
'ARS-USMARC-Parent-DQ888309-rs29013741': {'803', '804'},
'ARS-USMARC-Parent-DQ888310-rs29012422': {'805', '806'},
'ARS-USMARC-Parent-DQ888311-rs29017313': {'807', '808'},
'ARS-USMARC-Parent-DQ888313-no-rs': {'809', '810'},
'ARS-USMARC-Parent-DQ916057-rs29009979': {'811', '812'},
'ARS-USMARC-Parent-DQ916058-rs29016146': {'813', '814'},
'ARS-USMARC-Parent-DQ916059-rs29009907': {'815', '816'},
'ARS-USMARC-Parent-DQ984825-rs29012457': {'817', '818'},
'ARS-USMARC-Parent-DQ984826-rs29027559': {'819', '820'},
'ARS-USMARC-Parent-DQ984827-rs29012019': {'821', '822'},
'ARS-USMARC-Parent-DQ984828-rs29010004': {'823', '824'},
'ARS-USMARC-Parent-DQ990832-rs29015065': {'825', '826'},
'ARS-USMARC-Parent-DQ990833-rs29010147': {'827', '828'},
'ARS-USMARC-Parent-DQ990834-rs29013727': {'829', '830'},
'ARS-USMARC-Parent-DQ995976-no-rs': {'831', '832'},
'ARS-USMARC-Parent-DQ995977-rs29020834': {'833', '834'},
'ARS-USMARC-Parent-EF026084-rs29025380': {'835', '836'},
'ARS-USMARC-Parent-EF026085-rs29021607': {'837', '838'},
'ARS-USMARC-Parent-EF026086-rs29013660': {'839', '840'},
'ARS-USMARC-Parent-EF026087-rs29011643': {'841', '842'},
'ARS-USMARC-Parent-EF028073-rs29014953': {'843', '844'},
'ARS-USMARC-Parent-EF034080-rs29024749': {'845', '846'},
'ARS-USMARC-Parent-EF034081-rs29009668': {'847', '848'},
'ARS-USMARC-Parent-EF034082-rs29013532': {'849', '850'},
'ARS-USMARC-Parent-EF034083-rs29018286': {'851', '852'},
'ARS-USMARC-Parent-EF034084-rs29016185': {'853', '854'},
'ARS-USMARC-Parent-EF034085-rs29025677': {'855', '856'},
'ARS-USMARC-Parent-EF034086-no-rs': {'857', '858'},
'ARS-USMARC-Parent-EF042090-no-rs': {'859', '860'},
'ARS-USMARC-Parent-EF042091-rs29014974': {'861', '862'},
'ARS-USMARC-Parent-EF089234-rs29020870': {'863', '864'},
'ARS-USMARC-Parent-EF093509-rs29015170': {'865', '866'},
'ARS-USMARC-Parent-EF093511-rs29012316': {'867', '868'},
'ARS-USMARC-Parent-EF093512-rs29013546': {'869', '870'},
'ARS-USMARC-Parent-EF141102-rs29015783': {'871', '872'},
'ARS-USMARC-Parent-EF150946-rs29023666': {'873', '874'},
'ARS-USMARC-Parent-EF164803-rs29011141': {'875', '876'},
'BTA-100337-no-rs': {'877', '878'},
'BTA-101724-no-rs': {'879', '880'},
'BTA-102818-no-rs': {'881', '882'},
'BTA-104712-no-rs': {'883', '884'},
'BTA-112867-no-rs': {'885', '886'},
'BTA-113124-no-rs': {'887', '888'},
'BTA-113829-no-rs': {'889', '890'},
'BTA-114831-no-rs': {'891', '892'},
'BTA-11701-rs29017459': {'893', '894'},
'BTA-118486-no-rs': {'895', '896'},
'BTA-119562-no-rs': {'897', '898'},
'BTA-122684-no-rs': {'899', '900'},
'BTA-14388-rs29023151': {'901', '902'},
'BTA-15593-no-rs': {'903', '904'},
'BTA-17762-no-rs': {'905', '906'},
'BTA-19852-no-rs': {'907', '908'},
'BTA-20841-no-rs': {'909', '910'},
'BTA-21671-no-rs': {'911', '912'},
'BTA-23353-no-rs': {'913', '914'},
'BTA-23892-no-rs': {'915', '916'},
'BTA-25146-no-rs': {'917', '918'},
'BTA-27953-no-rs': {'919', '920'},
'BTA-30857-no-rs': {'921', '922'},
'BTA-31397-no-rs': {'923', '924'},
'BTA-35637-no-rs': {'925', '926'},
'BTA-37062-no-rs': {'927', '928'},
'BTA-37834-no-rs': {'929', '930'},
'BTA-41091-no-rs': {'931', '932'},
'BTA-42041-no-rs': {'933', '934'},
'BTA-46150-no-rs': {'935', '936'},
'BTA-48320-no-rs': {'937', '938'},
'BTA-49375-no-rs': {'939', '940'},
'BTA-52296-no-rs': {'941', '942'},
'BTA-53892-no-rs': {'943', '944'},
'BTA-57495-no-rs': {'945', '946'},
'BTA-57610-no-rs': {'947', '948'},
'BTA-57925-no-rs': {'949', '950'},
'BTA-58073-no-rs': {'951', '952'},
'BTA-59053-no-rs': {'953', '954'},
'BTA-61758-no-rs': {'955', '956'},
'BTA-65013-no-rs': {'957', '958'},
'BTA-68995-no-rs': {'959', '960'},
'BTA-70433-no-rs': {'961', '962'},
'BTA-73768-no-rs': {'963', '964'},
'BTA-73797-no-rs': {'965', '966'},
'BTA-75252-no-rs': {'967', '968'},
'BTA-78066-no-rs': {'969', '970'},
'BTA-80863-no-rs': {'971', '972'},
'BTA-83844-no-rs': {'973', '974'},
'BTA-85566-no-rs': {'975', '976'},
'BTA-86052-no-rs': {'977', '978'},
'BTA-89424-no-rs': {'979', '980'},
'BTA-90783-no-rs': {'981', '982'},
'BTA-92021-no-rs': {'983', '984'},
'BTA-96274-no-rs': {'985', '986'},
'BTA-99659-no-rs': {'987', '988'},
'BTB-00012128': {'989', '990'},
'BTB-00052125': {'991', '992'},
'BTB-00061242': {'993', '994'},
'BTB-00076466': {'995', '996'},
'BTB-00079213': {'997', '998'},
'BTB-00082871': {'1000', '999'},
'BTB-00089384': {'1001', '1002'},
'BTB-00095699': {'1003', '1004'},
'BTB-00098773': {'1005', '1006'},
'BTB-00169573': {'1007', '1008'},
'BTB-00184322': {'1009', '1010'},
'BTB-00188171': {'1011', '1012'},
'BTB-00242750': {'1013', '1014'},
'BTB-00285653': {'1015', '1016'},
'BTB-00311926': {'1017', '1018'},
'BTB-00312091': {'1019', '1020'},
'BTB-00320041': {'1021', '1022'},
'BTB-00327755': {'1023', '1024'},
'BTB-00386041': {'1025', '1026'},
'BTB-00393790': {'1027', '1028'},
'BTB-00394801': {'1029', '1030'},
'BTB-00397198': {'1031', '1032'},
'BTB-00414664': {'1033', '1034'},
'BTB-00420215': {'1035', '1036'},
'BTB-00431734': {'1037', '1038'},
'BTB-00436535': {'1039', '1040'},
'BTB-00458773': {'1041', '1042'},
'BTB-00468476': {'1043', '1044'},
'BTB-00474688': {'1045', '1046'},
'BTB-00492128': {'1047', '1048'},
'BTB-00500979': {'1049', '1050'},
'BTB-00563896': {'1051', '1052'},
'BTB-00566358': {'1053', '1054'},
'BTB-00607392': {'1055', '1056'},
'BTB-00637941': {'1057', '1058'},
'BTB-00663735': {'1059', '1060'},
'BTB-00736933': {'1061', '1062'},
'BTB-00806621': {'1063', '1064'},
'BTB-00818821': {'1065', '1066'},
'BTB-00885200': {'1067', '1068'},
'BTB-00888771': {'1069', '1070'},
'BTB-00923512': {'1071', '1072'},
'BTB-00939306': {'1073', '1074'},
'BTB-00955215': {'1075', '1076'},
'BTB-00959009': {'1077', '1078'},
'BTB-00980953': {'1079', '1080'},
'BTB-01007411': {'1081', '1082'},
'BTB-01007746': {'1083', '1084'},
'BTB-01012231': {'1085', '1086'},
'BTB-01033227': {'1087', '1088'},
'BTB-01036181': {'1089', '1090'},
'BTB-01057979': {'1091', '1092'},
'BTB-01077379': {'1093', '1094'},
'BTB-01086841': {'1095', '1096'},
'BTB-01097609': {'1097', '1098'},
'BTB-01109703': {'1099', '1100'},
'BTB-01124378': {'1101', '1102'},
'BTB-01130079': {'1103', '1104'},
'BTB-01141508': {'1105', '1106'},
'BTB-01141770': {'1107', '1108'},
'BTB-01146838': {'1109', '1110'},
'BTB-01172317': {'1111', '1112'},
'BTB-01226542': {'1113', '1114'},
'BTB-01274700': {'1115', '1116'},
'BTB-01274755': {'1117', '1118'},
'BTB-01276763': {'1119', '1120'},
'BTB-01285245': {'1121', '1122'},
'BTB-01287574': {'1123', '1124'},
'BTB-01292634': {'1125', '1126'},
'BTB-01303828': {'1127', '1128'},
'BTB-01304704': {'1129', '1130'},
'BTB-01371672': {'1131', '1132'},
'BTB-01375460': {'1133', '1134'},
'BTB-01397485': {'1135', '1136'},
'BTB-01416427': {'1137', '1138'},
'BTB-01426876': {'1139', '1140'},
'BTB-01450068': {'1141', '1142'},
'BTB-01453354': {'1143', '1144'},
'BTB-01465034': {'1145', '1146'},
'BTB-01478115': {'1147', '1148'},
'BTB-01495784': {'1149', '1150'},
'BTB-01530236': {'1151', '1152'},
'BTB-01559194': {'1153', '1154'},
'BTB-01575668': {'1155', '1156'},
'BTB-01626709': {'1157', '1158'},
'BTB-01640085': {'1159', '1160'},
'BTB-01656080': {'1161', '1162'},
'BTB-01680332': {'1163', '1164'},
'BTB-01728863': {'1165', '1166'},
'BTB-01735218': {'1167', '1168'},
'BTB-01753605': {'1169', '1170'},
'BTB-01796032': {'1171', '1172'},
'BTB-01817680': {'1173', '1174'},
'BTB-01890990': {'1175', '1176'},
'BTB-01902778': {'1177', '1178'},
'BTB-01904985': {'1179', '1180'},
'BTB-01920914': {'1181', '1182'},
'BTB-01944037': {'1183', '1184'},
'BTB-01980499': {'1185', '1186'},
'BTB-01988933': {'1187', '1188'},
'BTB-02006134': {'1189', '1190'},
'Hapmap22963-BTA-91941': {'1191', '1192'},
'Hapmap23219-BTA-154903': {'1193', '1194'},
'Hapmap23726-BTC-051363': {'1195', '1196'},
'Hapmap23870-BTA-128381': {'1197', '1198'},
'Hapmap24215-BTA-163266': {'1199', '1200'},
'Hapmap24224-BTA-93350': {'1201', '1202'},
'Hapmap24414-BTC-073009': {'1203', '1204'},
'Hapmap25382-BTC-000577': {'1205', '1206'},
'Hapmap25708-BTC-043671': {'1207', '1208'},
'Hapmap26512-BTA-52638': {'1209', '1210'},
'Hapmap26649-BTA-128159': {'1211', '1212'},
'Hapmap26681-BTA-143845': {'1213', '1214'},
'Hapmap26728-BTA-153579': {'1215', '1216'},
'Hapmap26904-BTC-061747': {'1217', '1218'},
'Hapmap26973-BTA-147422': {'1219', '1220'},
'Hapmap27089-BTC-045849': {'1221', '1222'},
'Hapmap27294-BTC-032117': {'1223', '1224'},
'Hapmap30171-BTC-063976': {'1225', '1226'},
'Hapmap30362-BTA-79378': {'1227', '1228'},
'Hapmap30746-BTC-070031': {'1229', '1230'},
'Hapmap31098-BTA-136127': {'1231', '1232'},
'Hapmap32161-BTA-55429': {'1233', '1234'},
'Hapmap32642-BTA-160002': {'1235', '1236'},
'Hapmap32799-BTA-132450': {'1237', '1238'},
'Hapmap33032-BTA-148803': {'1239', '1240'},
'Hapmap33049-BTA-153946': {'1241', '1242'},
'Hapmap33373-BTA-148341': {'1243', '1244'},
'Hapmap33479-BTA-137049': {'1245', '1246'},
'Hapmap33504-BTA-154919': {'1247', '1248'},
'Hapmap33849-BES3_Contig195_423': {'1249', '1250'},
'Hapmap34413-BES11_Contig509_1565': {'1251', '1252'},
'Hapmap34424-BES10_Contig566_926': {'1253', '1254'},
'Hapmap34713-BES10_Contig634_914': {'1255', '1256'},
'Hapmap35220-BES9_Contig365_495': {'1257', '1258'},
'Hapmap35535-SCAFFOLD86180_8791': {'1259', '1260'},
'Hapmap35744-SCAFFOLD60587_8279': {'1261', '1262'},
'Hapmap35988-SCAFFOLD166019_6545': {'1263', '1264'},
'Hapmap36588-SCAFFOLD90561_9460': {'1265', '1266'},
'Hapmap38147-BTA-25404': {'1267', '1268'},
'Hapmap38187-BTA-105082': {'1269', '1270'},
'Hapmap38474-BTA-122598': {'1271', '1272'},
'Hapmap38899-BTA-52139': {'1273', '1274'},
'Hapmap39008-BTA-120199': {'1275', '1276'},
'Hapmap39103-BTA-111886': {'1277', '1278'},
'Hapmap39152-BTA-27712': {'1279', '1280'},
'Hapmap39375-BTA-112108': {'1281', '1282'},
'Hapmap39419-BTA-63237': {'1283', '1284'},
'Hapmap39425-BTA-70290': {'1285', '1286'},
'Hapmap39434-BTA-77640': {'1287', '1288'},
'Hapmap39461-BTA-109898': {'1289', '1290'},
'Hapmap39562-BTA-40351': {'1291', '1292'},
'Hapmap39690-BTA-82650': {'1293', '1294'},
'Hapmap39712-BTA-109811': {'1295', '1296'},
'Hapmap40148-BTA-92999': {'1297', '1298'},
'Hapmap40170-BTA-20573': {'1299', '1300'},
'Hapmap40508-BTA-21513': {'1301', '1302'},
'Hapmap40581-BTA-73771': {'1303', '1304'},
'Hapmap40696-BTA-24595': {'1305', '1306'},
'Hapmap40729-BTA-40319': {'1307', '1308'},
'Hapmap40735-BTA-43736': {'1309', '1310'},
'Hapmap41037-BTA-61760': {'1311', '1312'},
'Hapmap41097-BTA-85971': {'1313', '1314'},
'Hapmap41591-BTA-59790': {'1315', '1316'},
'Hapmap42104-BTA-121232': {'1317', '1318'},
'Hapmap42556-BTA-45815': {'1319', '1320'},
'Hapmap42588-BTA-56723': {'1321', '1322'},
'Hapmap42596-BTA-58793': {'1323', '1324'},
'Hapmap42605-BTA-61330': {'1325', '1326'},
'Hapmap42648-BTA-71195': {'1327', '1328'},
'Hapmap42774-BTA-104770': {'1329', '1330'},
'Hapmap43057-BTA-80741': {'1331', '1332'},
'Hapmap43073-BTA-83925': {'1333', '1334'},
'Hapmap43142-BTA-107561': {'1335', '1336'},
'Hapmap43297-BTA-57357': {'1337', '1338'},
'Hapmap43792-BTA-122725': {'1339', '1340'},
'Hapmap43815-BTA-23051': {'1341', '1342'},
'Hapmap43869-BTA-49512': {'1343', '1344'},
'Hapmap43906-BTA-66658': {'1345', '1346'},
'Hapmap43908-BTA-67758': {'1347', '1348'},
'Hapmap43953-BTA-83292': {'1349', '1350'},
'Hapmap44181-BTA-104938': {'1351', '1352'},
'Hapmap44357-BTA-24055': {'1353', '1354'},
'Hapmap44395-BTA-64864': {'1355', '1356'},
'Hapmap44400-BTA-72792': {'1357', '1358'},
'Hapmap44545-BTA-47658': {'1359', '1360'},
'Hapmap45337-BTA-109277': {'1361', '1362'},
'Hapmap45569-BTA-55944': {'1363', '1364'},
'Hapmap46411-BTA-15820': {'1365', '1366'},
'Hapmap46550-BTA-103548': {'1367', '1368'},
'Hapmap46653-BTA-47447': {'1369', '1370'},
'Hapmap46717-BTA-69803': {'1371', '1372'},
'Hapmap46758-BTA-108921': {'1373', '1374'},
'Hapmap46854-BTA-64068': {'1375', '1376'},
'Hapmap46894-BTA-89312': {'1377', '1378'},
'Hapmap47167-BTA-107276': {'1379', '1380'},
'Hapmap47185-BTA-114173': {'1381', '1382'},
'Hapmap47281-BTA-40051': {'1383', '1384'},
'Hapmap47624-BTA-44484': {'1385', '1386'},
'Hapmap47669-BTA-59022': {'1387', '1388'},
'Hapmap47849-BTA-117910': {'1389', '1390'},
'Hapmap47975-BTA-50226': {'1391', '1392'},
'Hapmap48018-BTA-60331': {'1393', '1394'},
'Hapmap48024-BTA-62291': {'1395', '1396'},
'Hapmap48222-BTA-122240': {'1397', '1398'},
'Hapmap48539-BTA-97781': {'1399', '1400'},
'Hapmap48543-BTA-98093': {'1401', '1402'},
'Hapmap48697-BTA-25943': {'1403', '1404'},
'Hapmap48717-BTA-31633': {'1405', '1406'},
'Hapmap49049-BTA-119736': {'1407', '1408'},
'Hapmap49061-BTA-122268': {'1409', '1410'},
'Hapmap49062-BTA-122353': {'1411', '1412'},
'Hapmap49083-BTA-21452': {'1413', '1414'},
'Hapmap49219-BTA-55225': {'1415', '1416'},
'Hapmap49333-BTA-82773': {'1417', '1418'},
'Hapmap49430-BTA-107417': {'1419', '1420'},
'Hapmap49452-BTA-112834': {'1421', '1422'},
'Hapmap49925-BTA-24427': {'1423', '1424'},
'Hapmap50033-BTA-56544': {'1425', '1426'},
'Hapmap50047-BTA-59192': {'1427', '1428'},
'Hapmap50048-BTA-59263': {'1429', '1430'},
'Hapmap50109-BTA-79983': {'1431', '1432'},
'Hapmap50113-BTA-80661': {'1433', '1434'},
'Hapmap50266-BTA-13664': {'1435', '1436'},
'Hapmap50486-BTA-83799': {'1437', '1438'},
'Hapmap50598-BTA-122724': {'1439', '1440'},
'Hapmap50837-BTA-98392': {'1441', '1442'},
'Hapmap51009-BTA-62629': {'1443', '1444'},
'Hapmap51227-BTA-41809': {'1445', '1446'},
'Hapmap51304-BTA-74559': {'1447', '1448'},
'Hapmap51305-BTA-74560': {'1449', '1450'},
'Hapmap51412-BTA-16926': {'1451', '1452'},
'Hapmap51418-BTA-20281': {'1453', '1454'},
'Hapmap51522-BTA-93536': {'1455', '1456'},
'Hapmap51527-BTA-97415': {'1457', '1458'},
'Hapmap51569-BTA-17751': {'1459', '1460'},
'Hapmap51598-BTA-47943': {'1461', '1462'},
'Hapmap51908-BTA-63031': {'1463', '1464'},
'Hapmap52022-BTA-40336': {'1465', '1466'},
'Hapmap52072-rs29018920': {'1467', '1468'},
'Hapmap52240-rs29013844': {'1469', '1470'},
'Hapmap52308-rs29009652': {'1471', '1472'},
'Hapmap52438-rs29009828': {'1473', '1474'},
'Hapmap52559-rs29015773': {'1475', '1476'},
'Hapmap52627-rs29014567': {'1477', '1478'},
'Hapmap53212-rs29015272': {'1479', '1480'},
'Hapmap53277-rs29025667': {'1481', '1482'},
'Hapmap53428-rs29019634': {'1483', '1484'},
'Hapmap54020-rs29023153': {'1485', '1486'},
'Hapmap54131-rs29019697': {'1487', '1488'},
'Hapmap54285-ss46526570': {'1489', '1490'},
'Hapmap54313-rs29012632': {'1491', '1492'},
'Hapmap54320-rs29012920': {'1493', '1494'},
'Hapmap54432-rs29022442': {'1495', '1496'},
'Hapmap54459-rs29016219': {'1497', '1498'},
'Hapmap54547-rs29012198': {'1499', '1500'},
'Hapmap54638-rs29022155': {'1501', '1502'},
'Hapmap54766-rs29012658': {'1503', '1504'},
'Hapmap54847-rs29022340': {'1505', '1506'},
'Hapmap55441-rs29010990': {'1507', '1508'},
'Hapmap55445-rs29014484': {'1509', '1510'},
'Hapmap56777-ss46526772': {'1511', '1512'},
'Hapmap56893-rs29013584': {'1513', '1514'},
'Hapmap57137-rs29011525': {'1515', '1516'},
'Hapmap57299-ss46527085': {'1517', '1518'},
'Hapmap57430-rs29016504': {'1519', '1520'},
'Hapmap57586-rs29015232': {'1521', '1522'},
'Hapmap57648-rs29022376': {'1523', '1524'},
'Hapmap57725-rs29022829': {'1525', '1526'},
'Hapmap58102-rs29011445': {'1527', '1528'},
'Hapmap58509-rs29024139': {'1529', '1530'},
'Hapmap58587-ss46526997': {'1531', '1532'},
'Hapmap58613-rs29012081': {'1533', '1534'},
'Hapmap58894-rs29013940': {'1535', '1536'},
'Hapmap58990-rs29020862': {'1537', '1538'},
'Hapmap59058-rs29016195': {'1539', '1540'},
'Hapmap59303-rs29022511': {'1541', '1542'},
'Hapmap59420-ss46527113': {'1543', '1544'},
'Hapmap59828-rs29027014': {'1545', '1546'},
'Hapmap59876-rs29018046': {'1547', '1548'},
'Hapmap60017-rs29023471': {'1549', '1550'},
'Hapmap60125-rs29018949': {'1551', '1552'},
'Hapmap60420-rs29010136': {'1553', '1554'},
'Hapmap60671-rs29023031': {'1555', '1556'},
'Hapmap60688-rs29014090': {'1557', '1558'},
'Hapmap60813-rs29010178': {'1559', '1560'},
'INRA-638': {'1561', '1562'},
'UA-IFASA-1517': {'1563', '1564'},
'UA-IFASA-4359': {'1565', '1566'},
'UA-IFASA-4515': {'1567', '1568'},
'UA-IFASA-4619': {'1569', '1570'},
'UA-IFASA-4785': {'1571', '1572'},
'UA-IFASA-4798': {'1573', '1574'},
'UA-IFASA-4904': {'1575', '1576'},
'UA-IFASA-5034': {'1577', '1578'},
'UA-IFASA-5633': {'1579', '1580'},
'UA-IFASA-6154': {'1581', '1582'},
'UA-IFASA-6532': {'1583', '1584'},
'UA-IFASA-7236': {'1585', '1586'},
'UA-IFASA-7925': {'1587', '1588'},
'UA-IFASA-7987': {'1589', '1590'},
'UA-IFASA-8224': {'1591', '1592'},
'UA-IFASA-8658': {'1593', '1594'},
'UA-IFASA-8833': {'1595', '1596'},
'UA-IFASA-8930': {'1597', '1598'},
'UA-IFASA-9039': {'1599', '1600'}}

class genZipPackage:
    def __init__(self, zipDatoteka):
        self.zipFile=zipfile.ZipFile(zipDatoteka)
        self.zipname=zipDatoteka
        self.name=zipDatoteka.strip(".zip")
        self.sernum=zipDatoteka.strip(".zip").strip('Matija_Rigler_')
        self.genodate=zipDatoteka.strip(".zip").strip("Matija_Rigler_")[-9:].strip("_").strip('/')
        self.infiles=self.zipFile.namelist()
        self.finalreportname=filter(lambda x: x.endswith ("FinalReport.zip"), self.infiles)
            
    def unzip(self):
        self.zipFile.extractall()
    
    def extractFinalReport(self): #extracts all FinalReports and extracts them / if there is no FinalReport, it returns notice
        if self.finalreportname:
            for i in self.finalreportname:
                self.zipFile.extract(i) 
                zipfile.ZipFile(i).extractall()
                os.remove(os.getcwd()+'/'+i)
        else:
            return 'No FinalReport infile in ' + self.name
        
    def extractSampleMap(self):
        if 'Sample_Map.zip' in self.infiles:
            self.zipFile.extract("Sample_Map.zip")           
            zipfile.ZipFile("Sample_Map.zip").extractall()
            os.remove('Sample_Map.zip')
            shutil.move('Sample_Map.txt', self.name+'_Sample_Map.txt')
        else:
            return 'No Sample_Map.zip in {} infiles'.format(self.name)
        
    def extractSampleNames(self):
        if 'Sample_Map.zip' in self.infiles: 
            self.extractSampleMap()
            #os.system('sed -i "s|SI  |SI|g" ' + self.name+'_Sample_Map.txt') #remove double spacing
            #os.system('sed -i "s|SI |SI|g" ' + self.name+'_Sample_Map.txt') #remove space in sample IDs
            sampleTable = pd.read_table(self.name+'_Sample_Map.txt')
            return list(sampleTable['Name'])
        else:
            return 'No Sample_Map.zip in {} infiles'.format(self.name)

    
    def extractErrorNames(self): #this creates a list of tuples - errorID, correctedID
        if 'Sample_Map.zip' in self.infiles:
            self.extractSampleMap()
            #os.system('sed -i "s|SI  |SI|g" ' + self.name+'_Sample_Map.txt') #remove double spacing
            #os.system('sed -i "s|SI |SI|g" ' + self.name+'_Sample_Map.txt') #remove space in sample IDs
            sampleTable = pd.read_table(self.name+'_Sample_Map.txt')
            names=sampleTable['Name']
            errornames=[]
            for i in names:
                try:
                    errornames.append((int(i), 'SI'+str(i)))
                except:
                    if i.find("  ") != -1:
                        a = i.replace("  ", " ")
                        if a.find(" ") != -1:
                            errornames.append((i,( "".join(a.split(" ")[0:2]))))
                    if i.find("  ") == -1 and i.find(" ") != -1:
                        errornames.append((i, "".join(i.split(" ")[:2])))
            return errornames
        else:
            return ('No Sample_Map.zip in {} infiles'.format(self.name))
                    
                  
    def extractSNPMap(self):
        if 'SNP_Map.zip' in self.infiles:
            self.zipFile.extract("SNP_Map.zip")            
            zipfile.ZipFile("SNP_Map.zip").extractall()
            os.remove('SNP_Map.zip')
            shutil.move('SNP_Map.txt', self.name+'_SNP_Map.txt')
        else:
            return 'No SNP_Map.zip in {} infiles'.format(self.name)

    def checkSubDir(self):
        if any(ime.endswith("/") for ime in self.infiles):
            return "Zip file contains multiple genotype packages, subdirectories will be created at file extraction"
            
    def zipSubDir(self, rmOriginalZip): #zips the 6 files and removes the directory
        if self.checkSubDir():
            self.subDirNames=filter(lambda x: x.endswith("/"), self.infiles)
            self.unzip()
            for i in self.subDirNames:
                shutil.make_archive(i, "zip",os.getcwd()+'/'+i)
                shutil.rmtree(os.getcwd()+'/'+i)
            
            if rmOriginalZip=='Y':
                os.remove(self.zipname)
                
                                     
                                                          
                                                                                                    
class pedFile:
    def __init__(self, pedDatoteka):
        self.pedname=pedDatoteka
        self.name=pedDatoteka.strip(".ped")
        self.pedContent=open(pedDatoteka).read().strip("\n").split("\n") #here don't read it in as panda table since it takes much longer
        self.pedContent=open(pedDatoteka).read().strip("\n").split("\n") #here don't read it in as panda table since it takes much longer
        self.samples=[line.split(" ")[1] for line in self.pedContent]
        self.mapContent=open(pedDatoteka.strip(".ped") + ".map").read().strip("\n").split("\n")
        try:
            self.snps=[line.split(" ")[1] for line in self.mapContent]
        except:
            self.snps=[line.split("\t")[1] for line in self.mapContent]
        try:
            self.chip=chips[(len(self.pedContent[0].split(" "))-6)/2]
        except: 
            self.chip = len(self.snps)
        
    
    def extractSNP(self, SNP):
        os.system("plink --file " + self.name + " --cow --extract-snp "+ SNP + " --recode --out " + SNP)
        
    def extractSNPList(self, SNPList, outName):
        if "SNPList.txt" in os.listdir(os.getcwd()):
            overwrite=raw_input("Existing SNPList.txt file in the current working directory.\
 Do you want to overwrite and proceed? \n Ped and Map files will also be overwritten. [Y/N] ")
            if overwrite=='N':
                return None
            if overwrite == 'Y':
                pass
        with open("SNPList.txt", "w") as f:
            for snp in SNPList:
                f.write(snp + "\n")
        os.system("plink --file " + self.name + " --cow --extract SNPList.txt --recode --out " + outName)

    def extractSNPTxt(self, SNPFile, outName):
        os.system("plink --file " + self.name + " --cow --extract " + SNPFile + " --recode --out " + outName)

    def extractSNPList_Binary(self, SNPList, outName):
        if "SNPList.txt" in os.listdir(os.getcwd()):
            overwrite=raw_input("Existing SNPList.txt file in the current working directory.\
 Do you want to overwrite and proceed? \n Ped and Map files will also be overwritten. [Y/N] ")
            if overwrite=='N':
                return None
            if overwrite == 'Y':
                pass
        with open("SNPList.txt", "w") as f:
            for snp in SNPList:
                f.write(snp + "\n")
        os.system("plink --file " + self.name + " --cow --extract SNPList.txt --make-bed --out " + outName)

    def extractTraitSNPs(self, Trait):
        SNPonChip = [x for x in TraitSNPs[Trait] if x in self.snps]
        if SNPonChip:
            print "\n"*3 +  "{0} SNPs found on the chip({1}): {2}".format(Trait, len(SNPonChip), SNPonChip)
            with open(Trait + ".txt", "w") as f:
                for snp in TraitSNPs[Trait]:
                    f.write(snp + "\n")
            os.system("plink --file " + self.name + " --cow --extract "+Trait+".txt --recode --out " + Trait)
        if not SNPonChip:
            return "No {} SNPs found on chip".format(Trait)
            
    def extractParentalSNPs(self, number):
        if number in (100,200,800):
            with open("ParentalSNP_" + str(number) + ".txt", "w") as f:
                for i in eval("SNP"+str(number)+"Sifrant_Dict"):
                    f.write(i + "\n")
            os.system("plink --file " + self.name + " --cow --extract ParentalSNP_800.txt --recode --out ParentalSNP800_"+self.sernum)
        else:
            return "Non-standard number of SNPs for parental verification"
  
    def individualSNPsDF(self, sampleID):
        position=[i for i,x in enumerate(self.samples) if x == sampleID][0]
        return  pd.DataFrame({self.pedContent[position].split(" ")[1] : self.pedContent[position].split(" ")[6:]})
        
    def individualSNPsList(self, sampleID):
        position=[i for i,x in enumerate(self.samples) if x == sampleID][0]
        return self.pedContent[position].split(" ")[6:]
        
class Concordance():
    def __init__(self, RefPed, CompPed):
        self.RefPed = pedFile(RefPed)
        self.CompPed = pedFile(CompPed)

    def concSNPs(self, SNPlist, outConcName):
        self.RefPed.extractSNPTxt(SNPlist, 'Ref_current')
        self.CompPed.extractSNPTxt(SNPlist, 'Comp_current')
        os.system(
            'plink --file Ref_current --merge Comp_current.ped Comp_current.map --merge-mode 7 --cow --out ' + outConcName)
        os.system('grep "for a concordance rate" ' + outConcName + '.log > '+outConcName+'.txt')

    def extractConc(self, ConcFile):
        return float(open(ConcFile, 'r').read().strip('.\n').split(" ")[-1])

class mapFile:
    def __init__(self, mapDatoteka):
        self.mapname=mapDatoteka
        self.name=mapDatoteka.strip(".map")
        self.sernum=mapDatoteka.strip(".map").strip("Matija_Rigler_")
        self.mapContent=pd.read_table(mapDatoteka, header=None, sep=" ")
        if len(self.mapContent.columns) != 4:
            self.mapContent=pd.read_table(mapDatoteka, header=None, sep="\t")
        self.mapContent.columns = ["chr", "snp", "start", "stop"]
        self.snps=self.mapContent["snp"]
        try:
            self.chip=chips[(len(self.mapContent))]       
        except:
            self.chip=len(self.snps)
        
    def chrSNPs(self, chromosome):
        printSNPs = raw_input("The number of SNPs on chromosome {0} is {1}. ".format(chromosome, len(self.mapContent[self.mapContent.chr==str(chromosome)]))+"Do you want to display all the SNPs? [Y/N] ")
        if printSNPs == 'Y':
            return self.mapContent[self.mapContent.chr==str(chromosome)]
        else:
            pass
            
            
    def posSNP(self, chromosome, pos):
        return self.mapContent[(self.mapContent.chr==str(chromosome)) & (self.mapContent.stop==pos)]
              
        
class mapFileMerged:
    def __init__(self, mapDatoteka):
        self.mapname=mapDatoteka
        self.name=mapDatoteka.strip(".map")
        self.sernum=mapDatoteka.strip(".map").strip("Matija_Rigler_")
        self.mapContent=pd.read_table(mapDatoteka, header=None, sep=" ", names=("chr", "snp", "start", "stop"))
        self.snps=self.mapContent["snp"]
        try:
            self.chip=chips[(len(self.mapContent))]       
        except:
            self.chip=len(self.snps)
        
    def chrSNPs(self, chromosome):
        printSNPs = raw_input("The number of SNPs on chromosome {0} is {1}. ".format(chromosome, len(self.mapContent[self.mapContent.chr==str(chromosome)]))+"Do you want to display all the SNPs? [Y/N] ")
        if printSNPs == 'Y':
            return self.mapContent[self.mapContent.chr==str(chromosome)]
        else:
            pass
            
            
    def posSNP(self, chromosome, pos):
        return self.mapContent[(self.mapContent.chr==str(chromosome)) & (self.mapContent.stop==pos)]


def mergeList(RefFile, List):
   with open("MergeLIST.txt", "w") as f:
       for i in List:
           f.write(i + "\n")
   os.system("plink --file " + RefFile + " --cow --merge-list MergeLIST.txt --recode --out MERGEDforImputation")

def mergeList_autosomes(RefFile, List):
   with open("MergeLIST.txt", "w") as f:
       for i in List:
           f.write(i + "\n")
   os.system("plink --file " + RefFile + " --cow --chr 1-29 --merge-list MergeLIST.txt --recode --out MERGEDforImputation")