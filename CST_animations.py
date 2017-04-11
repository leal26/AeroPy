#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np
import types

#create image with format (time,x,y)
image = np.random.rand(10,10,10)
image2 = np.random.rand(10,10,10)

#setup figure
fig = plt.figure()
ax1=fig.add_subplot(1,2,1)
ax2=fig.add_subplot(1,2,2)

#set up list of images for animation
ims=[]
for time in xrange(np.shape(image)[1]):
    im = ax1.imshow(image[time,:,:])
 #   im2 = ax2.imshow(image2[time,:,:])
    im2, = ax2.plot(image[0:time,5,5])
 #   def setvisible(self,vis):
 #       for c in self.collections: c.set_visible(vis)
 #    im2.set_visible = types.MethodType(setvisible,im2,None)
 #    im2.axes = plt.gca()

    ims.append([im, im2])

#run animation
ani = anim.ArtistAnimation(fig,ims, interval=50,blit=False)
plt.show()
