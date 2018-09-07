import SimpleITK as sitk
import os
os.chdir("/home/jana/Documents/PhD/ZajemInAnalizaSlike/Vaja3")

images = [img for img in os.listdir(os.getcwd()) if img.endswith(".nrrd")]
print(images)

for img in images:
    Img = sitk.ReadImage(img)
    sitk.WriteImage(Img, img.strip("nrrd") + "nii", True)



