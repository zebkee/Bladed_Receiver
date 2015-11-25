# Only does the raytracing renderer
# Import Zeb's functions and models
from functions.zebfn import *
from models.receiver import*
import time

# Define the receiver and parameters
absorp = 0.7
blades = 4
span = 1.
divs = 25.
E = 2.0

width = 1.
height = 1.
depth = 0.7
ratio = 0.5

a,b,c = true_bladed_receiver(absorp,blades,width,width,height,depth,ratio)

rays_per_hstat = 100
no_of_hstat = 218
min_rays = 1000000

raycount = 0
totalenergy = 0
hitraycount = 0

# Raytracing commences, first iteration shows the target receiver
scene = TowerScene(reclist = a, recobj = b, objlist = c, z_offset = 6.1)
scene.aim_field()
scene.trace(hstat_rays=1,render=True)
H,boundlist,extent,binarea = scene.z_histdata(50)
raycount += 1*no_of_hstat
energy,subhits = scene.z_totalenergy(absorp)

print("Iter 1")
print("Rays incident = "+str(raycount))
print("Hit rays = "+str(subhits))
print("Rays absorbed = "+str(energy))

n=2
while raycount < min_rays:
	print("Iter "+str(n))
	scene = TowerScene(reclist = a, recobj = b, objlist = c, z_offset = 6.1)
	scene.aim_field()
	scene.trace(hstat_rays=rays_per_hstat,render=False)
	H2, boundlist, extent, area = scene.z_histdata(50)
	H += H2
	raycount += rays_per_hstat*no_of_hstat
	energy,subhits = scene.z_totalenergy(absorp)
	hitraycount += subhits
	print("Rays incident = "+str(raycount))
	print("Hit rays = "+str(subhits))
	print("Rays absorbed = "+str(energy))
	n += 1

print("Total number of incoming rays is "+str(raycount))
print("Total number of on-target rays is "+str(subhits))
print("Total amount of rays absorbed is "+str(energy))

# State the spillage percentage
spill = 1.0*(raycount-subhits)/(raycount)
print("The percentage of spilt rays is "+str(100*spill)+" %")
# State the surface absorptivity
print("The surface absorptivity is "+str(100*absorp)+" %")
# Now show the effective absorptivity
#energy,subhits = scene.z_totalenergy(absorp)
effabs = energy/(subhits)
print("The effective absorptivity is "+str(100*effabs)+" %")
impratio = effabs/absorp
print("The improvement ratio is "+str(impratio))

# Now get a visual output
# Normalize the histogram energy weights
e_per_ray = 218.0/raycount #this is in Watts
#print("The energy per 1.0 unit of rays is "+str(e_per_ray)+" kW")
factor = e_per_ray/binarea
#print("The area of each bin is "+str(binarea)+" sq metres")
#print("Color bar of histogram gives units of kW/sqmetre")

# Generate the 2D histogram plot
H_norm = H*(factor)
img = plt.imshow(H_norm,extent=extent, interpolation='nearest')
plt.colorbar()
y = r_[0,height]
n = 1
for bound in boundlist:
	globals()['line%s' % n] = plt.plot(r_[bound,bound],y,color='k') #Shows plate boundaries
	n += 1
plt.xlim(0,boundlist[-1]) #bounds the graph
plt.show()

# Generate the 1D histogram plot
z_hist1d(H_norm,100)

#scene.z_eff()
#tic = time.clock()
#scene.z_box_approx(1.0,50)
#scene.z_box_approx_2()
#toc = time.clock()
#scene.z_renderV(E)
#print("time here")
#print(toc-tic)
