from __future__ import division
from __future__ import print_function
import numpy as np
import PIL
from PIL import Image
import matplotlib.pyplot as plt
import SimpleITK as sitk
import os
import pandas as pd


os.chdir("/home/jana/Documents/PhD/ZajemInAnalizaSlike/Vaja4")

ct = sitk.ReadImage("ct.nrrd", sitk.sitkUInt8)
xray = sitk.ReadImage("xray.nrrd", sitk.sitkUInt16)

ctNP = sitk.GetArrayFromImage(ct)
xrayNP = sitk.GetArrayFromImage(xray)
# plt.figure()
# plt.imshow(xrayNP)
# plt.show()
print(xrayNP)
print(ctNP.shape)

