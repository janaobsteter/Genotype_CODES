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


