import os
import scipy.misc
import shutil

for specie in ['Plant', 'Animal', 'Fungi']:
    os.chdir("/home/jana/Pictures/DNN/" + specie)
    number = 0
    for image in os.listdir(os.getcwd()):
        a = scipy.misc.imread(image)
        a = scipy.misc.imresize(a, (299, 299))
        scipy.misc.imsave("Image" + str(number) + ".jpg", a)
        number = number + 1



for specie in ['Plant', 'Animal', 'Fungi']:
    os.chdir("/home/jana/Pictures/DNN/" + specie)
    for number in range(1, 20):
        try:
            shutil.move("Image" + str(number) + ".jpg", species + str(number) + ".jpg")
        except:
            pass
