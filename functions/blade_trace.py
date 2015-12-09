from functions.tower_scene import *
from models.receiver import *
import matplotlib.pyplot as plt
import timeit

import copy

def generic_trace(index,rec_obj,surf_ls,crit_ls,absorp,rays_per_bund=10000,tot_rays=100000,bins=50,render_bool=True,T=[0.,0.,6.1,0.,0.,0.]):
	"""Raytraces a generic receiver, where rec_obj is an assembled object, surf_ls is a list of all surfaces and
	   crit_ls is a list of all surfaces that are to be plotted in the histogram
	   T is a list of transformations [dx,dy,dz,rx,ry,rz] to move or rotate the receiver"""
	heliocsv = "/home/zebedee/Tracer/Tracer/examples/sandia_hstat_coordinates.csv"
	start = timeit.default_timer()
	
	A = absorp

	#Creates copies of these so a to not affect the original setup
	rec_obj_copy = copy.deepcopy(rec_obj)
        surf_ls_copy = copy.deepcopy(surf_ls)
        crit_ls_copy = copy.deepcopy(crit_ls)

	rays_used = 0
	hitsT = 0
	energyT = 0

	rph = int(rays_per_bund/218)

	scene = TowerScene(rec_obj=rec_obj_copy,surf_ls = surf_ls_copy, crit_ls = crit_ls_copy, heliostat = heliocsv,\
		dx=T[0],dy=T[1],dz=T[2],rx=T[3],ry=T[4],rz=T[5])
	scene.aim_field()

	print("Start Tracing")
	
	scene.trace(rph=rph,render = render_bool)
	H1, boundlist, extent, area = scene.hist_comb(no_of_bins=bins)
	energy,hits = scene.energies(A)
	energyT = energy
	hitsT = hits
	binarea = area #The first assignment is enough
	#H2 = H1 #The first assignment is enough
	extent1 = extent
	boundlist1 = boundlist # First assignment is enough
	rays_used += rph*218
	del scene #Try to free up memory

	while rays_used < tot_rays:
		scene = TowerScene(rec_obj=rec_obj_copy,surf_ls = surf_ls_copy, crit_ls = crit_ls_copy, heliostat = heliocsv)
		scene.aim_field()
		scene.trace(rph=rph,render=False)
		rays_used += rph*218
	H1, boundlist, extent, area = scene.hist_comb(no_of_bins=bins)
	energy, hits = scene.energies(A)
	energyT = energy
	hitsT = hits
	H2 = H1
	del scene #Try to free up memory
	#Calculate the final results
	spill = 1.0*(rays_used-hitsT)/rays_used
	effabs = energyT/hitsT
	e_per_ray = 218.0/rays_used
	H_norm = H2*(e_per_ray/binarea)
	#print(H_norm)
	max_conc = N.amax(H_norm)
	stop = timeit.default_timer()
	time_used = stop-start
	print("Finished "+str(rays_used)+" rays")

	return [index,spill,effabs,max_conc,rays_used,time_used,H_norm,boundlist1,extent1]
	
def write_data(filename,lsofls):
	"""Writes data onto a given filename. lsofls refers to a list of list of results
	   such as [[0,1,2,3,4,5],[0,1,2,3,4,5]]. Normally the histogram and boundlist 
	   should be left out"""
	f = open(filename,'w')
	f.write("Index "+"Spill "+"EffAbs "+"MaxConc "+"TotRays "+"TotTime "+"\n")
	for line in lsofls:	#take the first 6 values only
		f.write(str(line[0])+" "+str(line[1])+" "+str(line[2])+" "+str(line[3])+\
		" "+str(line[4])+" "+str(line[5])+"\n")
	f.close

def show_hist(H_norm,boundlist,height,extent):
	"""Displays a histogram and the lines which separate each surface"""
	img = plt.imshow(H_norm,extent = extent,interpolation='nearest')
	plt.colorbar()

	# Display the vertical lines
	y = r_[0,height]
	n = 1
	for bound in boundlist:
		globals()['line%s' % n] = plt.plot(r_[bound,bound],y,color='k') #Shows plate boundaries
		n += 1
	plt.xlim(0,boundlist[-1]) #bounds the graph
	plt.show()

def save_hist(H_norm,boundlist,height,extent,index=1,filename="histogram",dpi=1000):
    	"""Saves a histogram and the lines which separate each surface"""
    	plt.rcParams['axes.linewidth'] = 0.5 #set the value globally
    	img = plt.imshow(H_norm,extent = extent,interpolation='nearest')
    	plt.colorbar()

    	# Display the vertical lines
    	y = r_[0,height]
    	n = 1
    	for bound in boundlist:
        	globals()['line%s' % n] = plt.plot(r_[bound,bound],y,linewidth=0.2,color='k') #Shows plate boundaries
        	n += 1
    	plt.xlim(0,boundlist[-1]) #bounds the graph
    	newstring = filename+str(index)+".png"
    	plt.xlabel("width (m)")
    	plt.ylabel("height (m)")
    	plt.title(newstring)
    	plt.savefig(open(newstring, 'w'),dpi=dpi)
    	plt.close('all')	
		
	
