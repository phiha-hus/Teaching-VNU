# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 11:04:59 2020

@author: phiha

Image compression using Singular Value Decomposition (svd)
"""

import os
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt


#path = 'Sleeping_baby.jpg'
path = 'Baby_AJ.jpg'

img = Image.open(path)
s = float(os.path.getsize(path))/1000
print("Size(dimension): ",img.size)
plt.title("Original Image (%0.2f Kb):" %s)
plt.imshow(img)

imggray = img.convert('LA')
imgmat = np.array( list(imggray.getdata(band = 0)), float)
imgmat.shape = (imggray.size[1], imggray.size[0])
imgmat = np.matrix(imgmat)
plt.figure()
plt.imshow(imgmat, cmap = 'gray')
plt.title("Image after converting it into the Grayscale pattern")
plt.show()


print("After compression: ")
U, S, Vt = np.linalg.svd(imgmat) #single value decomposition

print("dimension of A is: ",np.size(imgmat))

for i in range(50, 300, 50):
    cmpimg = np.matrix(U[:, :i]) * np.diag(S[:i]) * np.matrix(Vt[:i,:])
    plt.imshow(cmpimg, cmap = 'gray')
    title = " Image after =  %s" %i
    plt.title(title)
    plt.show()
    result = Image.fromarray((cmpimg ).astype(np.uint8))
    result.save('compressed after = %s.jpg'%i)
    

