# Only does the raytracing renderer
# Import Zeb's functions and models
from functions.zebfn import *
from models.receiver import*
import time

# Define the receiver
absorp = 0.7
blades = 4
span = 1.
divs = 25.
E = 2.0

width = 1.
height = 1.
depth = 0.7
ratio = 0.5
#ratio is ratio of blade thickness to spacing


#a,b = test_receiver(1.0,1.0,absorp)
a,b = true_bladed_receiver(absorp,blades,width,width,height,depth,ratio)
#a,b = thin_bladed_receiver(1.0,1.0, absorp, 2)
#a,b = cavity_receiver(1.0,1.0,absorp)



scene = TowerScene(reclist = a, recobj = b, z_offset = 6.1)
scene.aim_field()
scene.trace(hstat_rays=20,render=True)
#scene.z_hist()
H = scene.z_histall(50)

z_hist1d(H)

scene.z_eff()
tic = time.clock()
#scene.z_box_approx(1.0,50)
#scene.z_box_approx_2()
toc = time.clock()
#scene.z_renderV(E)
print("time here")
print(toc-tic)
