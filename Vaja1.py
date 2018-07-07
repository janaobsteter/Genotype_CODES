from __future__ import division
from __future__ import print_function
import numpy as np
import PIL
from PIL import Image
#import matplotlib.pyplot as plt
import SimpleITK as sitk
import os
import pandas as pd

import scipy.ndimage as ni
#import matplotlib.pyplot as plt
from scipy.interpolate import interpn


plt.style.use('ggplot')





zob = PIL.Image.open(Dir + "zob-microct.png")
#zob.show()
print(zob.size) #get the size
print(zob.mode)
zob.getbands()

zobNP = np.array(zob)
#print(zobNP.view())
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

zobFloat = sitk.Cast(zobS, sitk.sitkFloat32 )
#nelinearno filtriranje
zobSc = sitk.GradientAnisotropicDiffusion(zobFloat, timeStep=0.125,
                                  conductanceParameter=int(2.0),
                                  conductanceScalingUpdateInterval=int(1),
                                  numberOfIterations=int(20))

ZOB = sitk.Cast( zobSc, sitk.sitkUInt8 )
sitk.WriteImage(ZOB, Dir + "Zob_Anisotropic1.jpg")


########################################################################################
########################################################################################
#DODATNE NALOGE

#control picture area
#homogoenous intensity area
#zob_control = zobNP[:150, :150]
PIL.Image.fromarray(zob_control).convert("L").show()

#ustvari preglednico, ki bo držala rezultate filtriranja
filtering = pd.DataFrame(columns=["Method", "Parameters", "SD"])

#testiraj različne vrednosti parametrov za anizortopni filter
for timeStep in [0.025, 0.125, 0.25]:
    for condParam in [1.0, 1.5, 2.0, 2.5, 3.0]:
        for scaling in range(1, 5):
            for iter in [5, 10, 15, 20, 25]:
                zobSc = sitk.CurvatureAnisotropicDiffusion(zobFloat, timeStep=timeStep,
                                                  conductanceParameter=int(condParam),
                                                  conductanceScalingUpdateInterval=int(scaling),
                                                  numberOfIterations=int(iter))

                ZOB = sitk.Cast( zobSc, sitk.sitkUInt8)
                #zapiši parametre filtiranja in preglej kakovost filtriranja s standardnim odklonom intenzitet na homogenem obmocju
                filtering = filtering.append(pd.DataFrame({"Method": ["curvature"],
                                                           "Parameters": (str(timeStep) + "_" + str(condParam) + "_" + str(scaling) + "_" + str(iter)),
                                                           "SD": [np.std(sitk.GetArrayFromImage(ZOB)[:150, :150])]}), sort=False)

#print(filtering)


#BILATERAL
#testiraj različne vrednosti parametrov za bilateralen filter
for samples in [1, 10, 20, 30, 50]:
    for domainSigma in [9, 18, 27]:
        #print("Mean", np.mean(np.gradient(zobNP)))
        #print("Median", np.median(np.gradient(zobNP)))
        zobBil_mean = sitk.Bilateral(zobFloat, domainSigma,  np.mean(np.gradient(zobNP)), samples)
        zobBil_Median = sitk.Bilateral(zobFloat, domainSigma,  np.median(np.gradient(zobNP)), samples)

        zobSc = sitk.Cast(zobBil_mean, sitk.sitkUInt8)
        ZOB = sitk.Cast(zobSc, sitk.sitkUInt8)
        # zapiši parametre filtiranja in preglej kakovost filtriranja s standardnim odklonom intenzitet na homogenem obmocju
        #filtriranje s povprečjem
        filtering = filtering.append(pd.DataFrame({"Method": ["Bilateral"],
                                                   "Parameters": ("Mean_" + str(samples) + "_" + str(domainSigma)),
                                                   "SD": [np.std(sitk.GetArrayFromImage(ZOB)[:150, :150])]}), sort=False)

        zobSc = sitk.Cast(zobBil_Median, sitk.sitkUInt8)
        ZOB = sitk.Cast(zobSc, sitk.sitkUInt8)
        # zapiši parametre filtiranja in preglej kakovost filtriranja s standardnim odklonom intenzitet na homogenem obmocju
        # filtriranje z mediano
        filtering = filtering.append(pd.DataFrame({"Method": ["Bilateral"],
                                                   "Parameters": ("Median_" + str(samples)+ "_" + str(domainSigma)),
                                                   "SD": [np.std(sitk.GetArrayFromImage(ZOB)[:150, :150])]}), sort=False)

#zapisi preglednico z rezultati
filtering.to_csv(Dir + "Filtering_Results.csv")

#shrani najuspešnejše rezultate filtriranja - slike
zobBil_mean = sitk.Bilateral(zobFloat, 18,  np.mean(np.gradient(zobNP)), 1)
ZOB = sitk.Cast( zobBil_mean, sitk.sitkUInt8 )
sitk.WriteImage(ZOB, Dir + "Zob_Bilateral_Mean1.jpg")




#%% Naloga 2
# # def getCubicBSpline2DGrid(iImageSize, iStep, nxy=None):
# #     return oCPx, oCPy
#
# # test funkcije
# if __name__ == '__main__':
#     cbImage = getChessBoardImage((400,400), 50)
#     iStep = (80, 80)
#     oCPx, oCPy = getCubicBSpline2DGrid(cbImage.shape, iStep)
#
#     plt.close('all')
#     plt.figure()
#     plt.imshow(cbImage, cmap='gray')
#     plt.plot(oCPx, oCPy, marker='o', color='r', linewidth=1)
#     plt.plot(oCPx.transpose(), oCPy.transpose(), marker='o', color='r', linewidth=1)
#     plt.xlim([-150,650])
#     plt.ylim([650,-150])
#     plt.show()
#
#
# #%% Naloga 3
# # def getCubicBSpline2DDeformation(iImageSize, iCPx, iCPy, iStep):
# #     return oGx, oGy
#
#
# # test funkcije
# if __name__ == '__main__':
#     cbImage = getChessBoardImage((400,400), 50)
#
#     iStep = (80, 80)
#     oCPx, oCPy = getCubicBSpline2DGrid(cbImage.shape, iStep)
#     oCPx_swirl, oCPy_swirl = swirlControlPoints(oCPx, oCPy, a=2.0, b=100.0)
#
#     oGx, oGy = getCubicBSpline2DDeformation(cbImage.shape, oCPx_swirl, oCPy_swirl, iStep)
#
#     # tri osi v enem prikaznem oknu
#     f, (ax1, ax2, ax3) = plt.subplots(1, 3)
#     ax1.imshow(cbImage, cmap='gray')
#     ax1.plot(oCPx_swirl, oCPy_swirl, marker='o', color='r', linewidth=1)
#     ax1.plot(oCPx_swirl.transpose(), oCPy_swirl.transpose(), marker='o', color='r', linewidth=1)
#     ax2.imshow(oGx, cmap='jet')
#     ax3.imshow(oGy, cmap='jet')
#
#
#
# #%% Naloga 4
# # def deformImageBSpline2D(iImage, iCPx, iCPy, iStep, interp_method='linear', inverse=False):
# #     return oImage
#
# # test funkcije
# if __name__ == '__main__':
#     cbImage = getChessBoardImage((400,400), 50)
#
#     iStep = (80, 80)
#     oCPx, oCPy = getCubicBSpline2DGrid(cbImage.shape, iStep)
#     oCPx_swirl, oCPy_swirl = swirlControlPoints(oCPx, oCPy, a=2.0, b=100.0)
#     cbImageDeformed = deformImageBSpline2D(cbImage, oCPx_swirl, oCPy_swirl, iStep)
#
#     # dvoje/troje osi v enem prikaznem oknu
#     f, (ax1, ax2) = plt.subplots(1, 2)
# #    f, (ax1, ax2, ax3) = plt.subplots(1, 3)
#     ax1.imshow(cbImage, cmap='gray')
#     ax2.imshow(cbImageDeformed, cmap='gray')
# #    ax3.imshow(cbImageDeformed, cmap='gray')
# #    ax3.plot(oCPx_swirl, oCPy_swirl, marker='o', color='r', linewidth=1)
# #    ax3.plot(oCPx_swirl.transpose(), oCPy_swirl.transpose(), marker='o', color='r', linewidth=1)
#     plt.show()
#
# #%% Naloga 5
# # test funkcije
# if __name__ == '__main__':
#     # nalozi slike
#     fixed = sitk.ReadImage('mr-nonenhanced.png', sitk.ssitkFloat32)
#     moving = sitk.ReadImage('mr-enhanced.png', sitk.ssitkFloat32)
#
#     # inicializacija postopka
#     R = sitk.ImageRegistrationMethod()
#
#     # inicializacija preslikave z B-zlepki
#     bsplineGrid = 8
#     bTr = sitk.BSplineTransformInitializer(fixed, [bsplineGrid] * 2)
#     R.SetInitialTransform(bTr, inPlace=True)
#
#     # inicializacija mere podobnosti
#     R.SetMetricAsMattesMutualInformation(50)
#     R.SetMetricSamplingPercentage(0.10)
#     R.SetMetricSamplingStrategy(R.RANDOM)
#
#     # inicializacija optimizacije
#     R.SetOptimizerAsGradientDescentLineSearch(learningRate=5.0,
#         numberOfIterations=100,
#         convergenceMinimumValue=1e-5,
#         convergenceWindowSize=5)
#     R.SetOptimizerScalesFromPhysicalShift()
#
#     # zagon poravnave
#     outTx = R.Execute(fixed, moving)
#
#     # ustvarjanje izhodne slike
#     S = sitk.ResampleImageFilter()
#     S.SetReferenceImage(fixed)
#     S.SetInterpolator(sitk.ssitkLinear)
#     S.SetDefaultPixelValue(0)
#     S.SetTransform(outTx)
#     outImage = S.Execute(moving)
#
#     # shranjevanje izhodne slike
#     sitk.WriteImage(sitk.Cast(outImage, sitk.ssitkUInt8),
#         'mr-enhanced-registered.png', True)
#
#     # prikaz slik pred in po poravnavi
#     fixed = sitk.ReadImage('mr-nonenhanced.png', sitk.ssitkUInt8)
#     moving = sitk.ReadImage('mr-enhanced.png', sitk.ssitkUInt8)
#     registered = sitk.ReadImage('mr-enhanced-registered.png', sitk.ssitkUInt8)
#
#     fa = sitk.GetArrayFromImage(fixed)
#     ma = sitk.GetArrayFromImage(moving)
#     ra = sitk.GetArrayFromImage(registered)
#
#     # šest osi v enem prikaznem oknu
#     f, ax = plt.subplots(2, 3)
#     ax[0,0].imshow(fa, cmap='gray')
#     ax[0,0].set_title('Fixed')
#     ax[0,1].imshow(ma, cmap='gray')
#     ax[0,1].set_title('Moving')
#     ax[0,2].imshow(ra, cmap='gray')
#     ax[0,2].set_title('Registered')
#
#     fa = normalizeImage(fa.astype('float'))
#     ma = normalizeImage(ma.astype('float'))
#     ra = normalizeImage(ra.astype('float'))
#
#     dy, dx = np.min(np.array((fa.shape,ma.shape)),axis=0)
#     ax[1,1].imshow(np.abs(ma[:dy,:dx]-fa[:dy,:dx]), cmap='gray')
#     ax[1,1].set_title('abs(Moving - Fixed)')
#     ax[1,2].imshow(np.abs(ra-fa), cmap='gray')
#     ax[1,2].set_title('abs(Registered - Fixed)')
#
#     plt.show()



# print(zobNP.view())
# print(np.size(zobNP, axis=1))
# print(np.size(zobNP, axis=0))
# print(zobNP[:2, :2])
#
# plt.figure()
# #plt.imshow(zobNP)
# #plt.show()
#
# #get the size if the image
# x = np.size(zobNP, axis=0)
# y = np.size(zobNP, axis=1)
#
# #half the x and y size
# zob_small = zobNP[:int(x/2), :int(y/2)]
# plt.imshow(zob_small)
# print(int(x/2), int(y/2))
# #plt.show()
#
# #save image
# PIL.Image.fromarray(zob_small).save(Dir + "zobSmall.jpg")
# #PIL.Image.fromarray(zob_small).convert("RGB").show()
#
#
# #vaja 1.1 - simpleITK
#
# #zobS = sitk.GetImageFromArray(zobNP)
# zobS = sitk.ReadImage(Dir + "zob-microct.png")
# print("Size: ", zobS.GetSize())
# print("Tip: ", zobS.GetPixelIDValue())
# print("Tip: ", zobS.GetPixelIDTypeAsString())
# print("Korak: ", zobS.GetSpacing())
# print("Koordinate: ", zobS.GetOrigin())
# print("Smer: ", zobS.GetDirection())
#
# #ponastavi korak vzorčenja
# zobSa = zobS #da ne spreminjam originalne slike
# zobSa.SetSpacing([2.0, 2.0])
# print("Nov korak vzorčenja: ", zobSa.GetSpacing())
# #ponastavi izhodišče
# zobSa.SetOrigin([3.0, 0.5])
# print("Novo izhodišče: ", zobSa.GetOrigin())
# sitk.WriteImage(zobSa, Dir + "Zob.nrrd")
# sitk.WriteImage(zobSa, Dir + "Zob.nii.gz")
#
# #naloži nrrd .nazaj
# zobS1 = sitk.ReadImage(Dir + "Zob.nrrd")
# print("Izhodišče .nrrd: ", zobS1.GetOrigin())
# zobS1_array = sitk.GetArrayFromImage(zobS1)
# print(zobS1_array.view())
#
# #spremeni iz array-a nazaj v sliko - preveri, če se je ohranilo izhodišče
# zobB = sitk.GetImageFromArray(zobS1_array)
# print("Izhodišče array: ", zobB.GetOrigin())
# print("Korak vzorčenja array: ", zobB.GetSpacing())
#
#
# #obnova medicinskih slik
# #filtriranje
# zobS_filter = sitk.Mean(zobS)
# sitk.WriteImage(zobS_filter, Dir + "Zob_filterMean.png")
#
# zobS_filter = sitk.Median(zobS)
# sitk.WriteImage(zobS_filter, Dir + "Zob_filterMedian.png")
#
# zobFloat = sitk.Cast( zobS, sitk.sitkFloat32 )
# #nelinearno filtriranje
# zobSc = sitk.GradientAnisotropicDiffusion(zobFloat, timeStep=0.125,
#                                   conductanceParameter=int(2.0),
#                                   conductanceScalingUpdateInterval=int(1),
#                                   numberOfIterations=int(20))
#
# ZOB = sitk.Cast( zobSc, sitk.sitkUInt8 )
# sitk.WriteImage(ZOB, Dir + "Zob_Anisotropic1.jpg")
#
#
# #########################################################################################
# #########################################################################################
# #DODATNE NALOGE
#
# #control picture area
# #homogoenous intensity area
# #zob_control = zobNP[:150, :150]
# #PIL.Image.fromarray(zob_control).convert("L").show()
# # filtering = pd.DataFrame(columns=["Method", "Parameters", "SD"])
# #
# # for timeStep in [0.025, 0.125, 0.25]:
# #     for condParam in [1.0, 1.5, 2.0, 2.5, 3.0]:
# #         for scaling in range(1, 5):
# #             for iter in [5, 10, 15, 20, 25]:
# #                 zobSc = sitk.CurvatureAnisotropicDiffusion(zobFloat, timeStep=timeStep,
# #                                                   conductanceParameter=int(condParam),
# #                                                   conductanceScalingUpdateInterval=int(scaling),
# #                                                   numberOfIterations=int(iter))
# #
# #                 ZOB = sitk.Cast( zobSc, sitk.sitkUInt8 )
# #                 #sitk.WriteImage(ZOB, Dir + "Zob_Curvature.jpg")
# #                 #print(np.std(sitk.GetArrayFromImage(ZOB)[:150, :150]))
# #                 filtering = filtering.append(pd.DataFrame({"Method": ["curvature"],
# #                                                            "Parameters": (str(timeStep) + "_" + str(condParam) + "_" + str(scaling) + "_" + str(iter)),
# #                                                            "SD": [np.std(sitk.GetArrayFromImage(ZOB)[:150, :150])]}), sort=False)
# #
# # print(filtering)
# #
# #
# # #BILATERAL
# # """
# # bilpd = pd.dataframe(columns=["method", "domainsigma", "rangesigma","gauss", "sd"])
# # for domainsigma in [0.1 * x for x in range (1, 20)]:
# #     for rangesigma in [0.1 * x for x in range (1, 20)]:
# #         for gauss in range(1, 20):
# #             zobsc = sitk.bilateral(zobfloat, domainsigma, rangesigma, 10)
# #
# #             zob = sitk.cast( zobsc, sitk.sitkuint8 )
# #             sitk.writeimage(zob, dir + "zob_bilateral" + str(domainsigma) + "_" + str(rangesigma) + ".jpg")
# #             bilpd= bilpd.append(
# #                 pd.dataframe({"method": ["bilateral"], "domainsigma": domainsigma,
# #                               "rangesigma": rangesigma, "gauss": gauss, "sd": [np.std(sitk.getarrayfromimage(zob)[:150, :150])]}), sort=false)
# # """
# #
# # for samples in [1, 10, 20, 30, 50]:
# #     for domainSigma in [9, 18, 27]:
# #         #print("Mean", np.mean(np.gradient(zobNP)))
# #         #print("Median", np.median(np.gradient(zobNP)))
# #         zobBil_mean = sitk.Bilateral(zobFloat, domainSigma,  np.mean(np.gradient(zobNP)), samples)
# #         zobBil_Median = sitk.Bilateral(zobFloat, domainSigma,  np.median(np.gradient(zobNP)), samples)
# #
# #         zobSc = sitk.Cast(zobBil_mean, sitk.sitkUInt8)
# #         ZOB = sitk.Cast(zobSc, sitk.sitkUInt8)
# #         filtering = filtering.append(pd.DataFrame({"Method": ["Bilateral"],
# #                                                    "Parameters": ("Mean_" + str(samples) + "_" + str(domainSigma)),
# #                                                    "SD": [np.std(sitk.GetArrayFromImage(ZOB)[:150, :150])]}), sort=False)
# #
# #         zobSc = sitk.Cast(zobBil_Median, sitk.sitkUInt8)
# #         ZOB = sitk.Cast(zobSc, sitk.sitkUInt8)
# #         filtering = filtering.append(pd.DataFrame({"Method": ["Bilateral"],
# #                                                    "Parameters": ("Median_" + str(samples)+ "_" + str(domainSigma)),
# #                                                    "SD": [np.std(sitk.GetArrayFromImage(ZOB)[:150, :150])]}), sort=False)
#
# # print(filtering)
# # filtering.to_csv(Dir + "Filtering_Results.csv")
#
# zobBil_mean = sitk.Bilateral(zobFloat, 18,  np.mean(np.gradient(zobNP)), 1)
# ZOB = sitk.Cast( zobBil_mean, sitk.sitkUInt8 )
# sitk.WriteImage(ZOB, Dir + "Zob_Bilateral_Mean1.jpg")
#
# zobBil_mean = sitk.Bilateral(zobFloat, 18,  np.median(np.gradient(zobNP)), 1)
# ZOB = sitk.Cast( zobBil_mean, sitk.sitkUInt8 )
# sitk.WriteImage(ZOB, Dir + "Zob_Bilateral_Median1.jpg")
#
# zobBil_mean = sitk.Bilateral(zobFloat, 18,  np.mean(np.gradient(zobNP)), 10)
# ZOB = sitk.Cast( zobBil_mean, sitk.sitkUInt8 )
# sitk.WriteImage(ZOB, Dir + "Zob_Bilateral_Mean10.jpg")
#
# zobBil_mean = sitk.Bilateral(zobFloat, 18,  np.mean(np.gradient(zobNP)), 20)
# ZOB = sitk.Cast( zobBil_mean, sitk.sitkUInt8 )
# sitk.WriteImage(ZOB, Dir + "Zob_Bilateral_Mean20.jpg")
#
# zobSc = sitk.CurvatureAnisotropicDiffusion(zobFloat, timeStep=0.25,
#                                   conductanceParameter=int(3.0),
#                                   conductanceScalingUpdateInterval=int(4),
#                                   numberOfIterations=int(25))
#

ZOB = sitk.Cast( zobSc, sitk.sitkUInt8 )
sitk.WriteImage(ZOB, Dir + "Zob_Curvature_Best.jpg")


filtering = pd.DataFrame(columns=["Method", "MaxLevel", "MaxIter" "SD"])

misice = sitk.ReadImage(Dir + "misice-microscope.png")
mis = PIL.Image.open(Dir + "misice-microscope.png")
misH = mis.histogram()
plt.hist(misH, bins=50)
#plt.show()
misiceNP = sitk.GetArrayFromImage(misice)
plt.hist(np.diagonal(misiceNP), bins=50)
plt.savefig(Dir + "Hist_Original.png")


#testiraj različne parametre za N4 postopek filtriranja
for iter in [10, 50, 100]:
    for level in range(1, 5):
        misice = sitk.ReadImage(Dir + "misice-microscope.png")
#        misice_np = np.array(misice)
        misiceFloat = sitk.Cast(misice, sitk.sitkFloat32)
        corrector = sitk.N4BiasFieldCorrectionImageFilter()
        corrector.SetMaximumNumberOfIterations([iter] * level)
        corrector.SetMaskLabel(1)
        oImage = corrector.Execute(misiceFloat)
        MISICE = sitk.Cast(oImage, sitk.sitkUInt8)
        MISICE_np = sitk.GetArrayFromImage(MISICE)
        #print(MISICE_np)
        plt.clf()
        plt.hist(np.diagonal(MISICE_np), bins=50)
        #plt.show()
        plt.savefig(Dir + "Hist_" + str(iter) + "_" + str(level) + ".png")



        sitk.WriteImage(MISICE, Dir + "MisiceN4" + str(iter) + "_" + str(level) + ".jpg")


