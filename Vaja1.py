from __future__ import division
from __future__ import print_function
import numpy as np
import PIL
import matplotlib.pyplot as plt
import SimpleITK as sitk
import os
import pandas as pd



Dir = "/home/jana/Documents/PhD/ZajemInAnalizaSlike/Vaja1/"

zob = PIL.Image.open(Dir + "zob-microct.png")
zob.show()
print(zob.size) #get the size
print(zob.mode)
zob.getbands()

zobNP = np.array(zob)
print(zobNP.view())
print(np.size(zobNP, axis=1))
print(np.size(zobNP, axis=0))
print(zobNP[:2, :2])

plt.figure()
#plt.imshow(zobNP)
#plt.show()

#get the size if the image
x = np.size(zobNP, axis=0)
y = np.size(zobNP, axis=1)

#half the x and y size
zob_small = zobNP[:int(x/2), :int(y/2)]
plt.imshow(zob_small)
print(int(x/2), int(y/2))
#plt.show()

#save image
PIL.Image.fromarray(zob_small).save(Dir + "zobSmall.jpg")
#PIL.Image.fromarray(zob_small).convert("RGB").show()


#vaja 1.1 - simpleITK

#zobS = sitk.GetImageFromArray(zobNP)
zobS = sitk.ReadImage(Dir + "zob-microct.png")
print("Size: ", zobS.GetSize())
print("Tip: ", zobS.GetPixelIDValue())
print("Tip: ", zobS.GetPixelIDTypeAsString())
print("Korak: ", zobS.GetSpacing())
print("Koordinate: ", zobS.GetOrigin())
print("Smer: ", zobS.GetDirection())

#ponastavi korak vzorčenja
zobSa = zobS #da ne spreminjam originalne slike
zobSa.SetSpacing([2.0, 2.0])
print("Nov korak vzorčenja: ", zobSa.GetSpacing())
#ponastavi izhodišče
zobSa.SetOrigin([3.0, 0.5])
print("Novo izhodišče: ", zobSa.GetOrigin())
sitk.WriteImage(zobSa, Dir + "Zob.nrrd")
sitk.WriteImage(zobSa, Dir + "Zob.nii.gz")

#naloži nrrd .nazaj
zobS1 = sitk.ReadImage(Dir + "Zob.nrrd")
print("Izhodišče .nrrd: ", zobS1.GetOrigin())
zobS1_array = sitk.GetArrayFromImage(zobS1)
print(zobS1_array.view())

#spremeni iz array-a nazaj v sliko - preveri, če se je ohranilo izhodišče
zobB = sitk.GetImageFromArray(zobS1_array)
print("Izhodišče array: ", zobB.GetOrigin())
print("Korak vzorčenja array: ", zobB.GetSpacing())


#obnova medicinskih slik
#filtriranje
zobS_filter = sitk.Mean(zobS)
sitk.WriteImage(zobS_filter, Dir + "Zob_filterMean.png")

zobS_filter = sitk.Median(zobS)
sitk.WriteImage(zobS_filter, Dir + "Zob_filterMedian.png")

zobFloat = sitk.Cast( zobS, sitk.sitkFloat32 )
#nelinearno filtriranje
zobSc = sitk.GradientAnisotropicDiffusion(zobFloat, timeStep=0.125,
                                  conductanceParameter=int(2.0),
                                  conductanceScalingUpdateInterval=int(1),
                                  numberOfIterations=int(20))

ZOB = sitk.Cast( zobSc, sitk.sitkUInt8 )
sitk.WriteImage(ZOB, Dir + "Zob_Anisotropic1.jpg")


#########################################################################################
#########################################################################################
#DODATNE NALOGE

#control picture area
#homogoenous intensity area
#zob_control = zobNP[:150, :150]
#PIL.Image.fromarray(zob_control).convert("L").show()
# filtering = pd.DataFrame(columns=["Method", "Parameters", "SD"])
#
# for timeStep in [0.025, 0.125, 0.25]:
#     for condParam in [1.0, 1.5, 2.0, 2.5, 3.0]:
#         for scaling in range(1, 5):
#             for iter in [5, 10, 15, 20, 25]:
#                 zobSc = sitk.CurvatureAnisotropicDiffusion(zobFloat, timeStep=timeStep,
#                                                   conductanceParameter=int(condParam),
#                                                   conductanceScalingUpdateInterval=int(scaling),
#                                                   numberOfIterations=int(iter))
#
#                 ZOB = sitk.Cast( zobSc, sitk.sitkUInt8 )
#                 #sitk.WriteImage(ZOB, Dir + "Zob_Curvature.jpg")
#                 #print(np.std(sitk.GetArrayFromImage(ZOB)[:150, :150]))
#                 filtering = filtering.append(pd.DataFrame({"Method": ["curvature"],
#                                                            "Parameters": (str(timeStep) + "_" + str(condParam) + "_" + str(scaling) + "_" + str(iter)),
#                                                            "SD": [np.std(sitk.GetArrayFromImage(ZOB)[:150, :150])]}), sort=False)
#
# print(filtering)
#
#
# #BILATERAL
# """
# bilpd = pd.dataframe(columns=["method", "domainsigma", "rangesigma","gauss", "sd"])
# for domainsigma in [0.1 * x for x in range (1, 20)]:
#     for rangesigma in [0.1 * x for x in range (1, 20)]:
#         for gauss in range(1, 20):
#             zobsc = sitk.bilateral(zobfloat, domainsigma, rangesigma, 10)
#
#             zob = sitk.cast( zobsc, sitk.sitkuint8 )
#             sitk.writeimage(zob, dir + "zob_bilateral" + str(domainsigma) + "_" + str(rangesigma) + ".jpg")
#             bilpd= bilpd.append(
#                 pd.dataframe({"method": ["bilateral"], "domainsigma": domainsigma,
#                               "rangesigma": rangesigma, "gauss": gauss, "sd": [np.std(sitk.getarrayfromimage(zob)[:150, :150])]}), sort=false)
# """
#
# for samples in [1, 10, 20, 30, 50]:
#     for domainSigma in [9, 18, 27]:
#         #print("Mean", np.mean(np.gradient(zobNP)))
#         #print("Median", np.median(np.gradient(zobNP)))
#         zobBil_mean = sitk.Bilateral(zobFloat, domainSigma,  np.mean(np.gradient(zobNP)), samples)
#         zobBil_Median = sitk.Bilateral(zobFloat, domainSigma,  np.median(np.gradient(zobNP)), samples)
#
#         zobSc = sitk.Cast(zobBil_mean, sitk.sitkUInt8)
#         ZOB = sitk.Cast(zobSc, sitk.sitkUInt8)
#         filtering = filtering.append(pd.DataFrame({"Method": ["Bilateral"],
#                                                    "Parameters": ("Mean_" + str(samples) + "_" + str(domainSigma)),
#                                                    "SD": [np.std(sitk.GetArrayFromImage(ZOB)[:150, :150])]}), sort=False)
#
#         zobSc = sitk.Cast(zobBil_Median, sitk.sitkUInt8)
#         ZOB = sitk.Cast(zobSc, sitk.sitkUInt8)
#         filtering = filtering.append(pd.DataFrame({"Method": ["Bilateral"],
#                                                    "Parameters": ("Median_" + str(samples)+ "_" + str(domainSigma)),
#                                                    "SD": [np.std(sitk.GetArrayFromImage(ZOB)[:150, :150])]}), sort=False)

# print(filtering)
# filtering.to_csv(Dir + "Filtering_Results.csv")

zobBil_mean = sitk.Bilateral(zobFloat, 18,  np.mean(np.gradient(zobNP)), 1)
ZOB = sitk.Cast( zobBil_mean, sitk.sitkUInt8 )
sitk.WriteImage(ZOB, Dir + "Zob_Bilateral_Mean1.jpg")

zobBil_mean = sitk.Bilateral(zobFloat, 18,  np.median(np.gradient(zobNP)), 1)
ZOB = sitk.Cast( zobBil_mean, sitk.sitkUInt8 )
sitk.WriteImage(ZOB, Dir + "Zob_Bilateral_Median1.jpg")

zobBil_mean = sitk.Bilateral(zobFloat, 18,  np.mean(np.gradient(zobNP)), 10)
ZOB = sitk.Cast( zobBil_mean, sitk.sitkUInt8 )
sitk.WriteImage(ZOB, Dir + "Zob_Bilateral_Mean10.jpg")

zobBil_mean = sitk.Bilateral(zobFloat, 18,  np.mean(np.gradient(zobNP)), 20)
ZOB = sitk.Cast( zobBil_mean, sitk.sitkUInt8 )
sitk.WriteImage(ZOB, Dir + "Zob_Bilateral_Mean20.jpg")

zobSc = sitk.CurvatureAnisotropicDiffusion(zobFloat, timeStep=0.25,
                                  conductanceParameter=int(3.0),
                                  conductanceScalingUpdateInterval=int(4),
                                  numberOfIterations=int(25))

ZOB = sitk.Cast( zobSc, sitk.sitkUInt8 )
sitk.WriteImage(ZOB, Dir + "Zob_Curvature_Best.jpg")