from functions.tower_scene import *
from models.receiver import *
import matplotlib.pyplot as plt
import timeit

import copy

def stats(array):
	"""Returns the mean and standard deviation percentage in decimal form of an array"""
	mean = float(sum(array))/len(array)
	std = (sum((array-mean)**2)/len(array))**0.5
	return mean,std

def generic_trace(index,rec_obj,surf_ls,crit_ls,absorp,rays_per_bund=10000,tot_rays=100000,bins=50,render_bool=True,T=[0.,0.,6.1,0.,0.,0.],sun_az = 0.,sun_elev = 34.9):
	"""Raytraces a generic receiver, where rec_obj is an assembled object, surf_ls is a list of all surfaces and
	   crit_ls is a list of all surfaces that are to be plotted in the histogram.
	   T is a transformation list of [dx,dy,dz,rx,ry,rz] to move and rotate the receiver"""
	heliocsv = "/home/zebedee-u5277975/Tracer/Tracer/examples/sandia_hstat_coordinates.csv"
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

	scene = TowerScene(rec_obj=rec_obj_copy,surf_ls = surf_ls_copy, crit_ls = crit_ls_copy, heliostat = heliocsv,dx=T[0],dy=T[1],dz=T[2],rx=T[3],ry=T[4],rz=T[5],sun_az=sun_az,sun_elev=sun_elev)
	scene.aim_field()
	
	scene.trace(rph=rph,render = render_bool)
	H1, boundlist, extent, area = scene.hist_comb(no_of_bins=bins)
	print(N.sum(H1))
	energy,hits,heliohits = scene.energies()
	energyT = energy
	hitsT = hits
	heliohitsT = heliohits
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
		energy, hits, heliohits = scene.energies()
		heliohitsT += heliohits
	energyT = energy
	hitsT = hits
		
	H2 = H1
	del scene #Try to free up memory
	#Calculate the final results
	spill = 1.0*(heliohitsT-hitsT)/heliohitsT
	effabs = energyT/hitsT
	e_per_ray = 37*218.0/heliohitsT
	H_norm = H2*(e_per_ray/binarea)
	#print(H_norm)
	max_conc = N.amax(H_norm)
	stop = timeit.default_timer()
	time_used = stop-start

	return [index,spill,effabs,max_conc,rays_used,time_used,H_norm,boundlist1,extent1]

def precision_trace(index,rec_obj,surf_ls,crit_ls,rays_per_iter=100000,bins=100,precision=2.0,render_bool=False,T=[0.,0.,6.1,0.,0.,0.],sun_az = 0.,sun_elev = 34.9):
	#precision is in %
	heliocsv = "/home/zebedee-u5277975/Tracer/Tracer/examples/sandia_hstat_coordinates.csv"
	print("99.7% confidence interval set to cover +/- "+str(precision)+ "% of the maxconc value")
	
	start = timeit.default_timer()

	rph = 1+int(rays_per_iter/218)
	rays_per_run = rph*218 #This is the actual number of rays per iteration
	total_rays = 0
	first_iter = True
	maxconc_array = N.array([]) #an empty array of MaxConc Values
	spill_array = N.array([])
	effabs_array = N.array([])
	maxconc_halfint_pct = 100.0 #This is the standard deviation of the maxconc value
	threshold = 0
	maxconc_mean = 1
	heliohitsT = 0

	ite = 1
	while maxconc_halfint_pct > precision or total_rays < 1000000 or ite < 10:
		print("Iteration no. "+str(ite))
		rec_obj_copy = copy.deepcopy(rec_obj)
        	surf_ls_copy = copy.deepcopy(surf_ls)
        	crit_ls_copy = copy.deepcopy(crit_ls)
		scene = TowerScene(rec_obj=rec_obj_copy,surf_ls = surf_ls_copy, crit_ls = crit_ls_copy,\
		        heliostat=heliocsv,dx=T[0],dy=T[1],dz=T[2],rx=T[3],ry=T[4],rz=T[5],sun_az=sun_az,sun_elev=sun_elev)
		scene.aim_field()
		scene.trace(rph=rph,render = render_bool)
		render_bool = False
		energy,hits,heliohits = scene.energies()
		heliohitsT += heliohits
		H1, boundlist, extent, area = scene.hist_comb(no_of_bins=bins)
		del scene
		if first_iter:
			H = H1
			first_iter = False
		else:
			H += H1
		
		maxconc = N.amax(H1)*37.0*218/(heliohits*area)
		maxconc_array = N.append(maxconc_array,maxconc)
		spill = 1.0*(heliohits-hits)/heliohits
		spill_array = N.append(spill_array,spill)
		effabs = float(energy)/hits
		effabs_array = N.append(effabs_array,effabs)
		
		maxconc_mean,maxconc_std = stats(maxconc_array)
		maxconc_halfint_pct = (100*(2.7478*maxconc_std)/(ite**0.5))/maxconc_mean
		print("Half interval is now "+str(maxconc_halfint_pct)+"%")  
		print("Maxconc = "+str(maxconc_mean)+" +/- "+str(maxconc_std))  
		spill_mean,spill_std = stats(spill_array)
		print("Spill = "+str(spill_mean)+" +/- "+str(spill_std))
		effabs_mean,effabs_std = stats(effabs_array)
		print("EffAbs = "+str(effabs_mean)+" +/- "+str(effabs_std))
		total_rays += rays_per_run
		ite += 1
		print("Processed "+str(total_rays)+" rays.")
		#print("Current %std is "+str(100*maxconc_std)+"%")
		print(" ")
	print("Satisfactory precision reached")
	#Final results
	maxconc_halfint_pct = (100*(2.7478*maxconc_std)/(ite**0.5))/maxconc_mean
	spill_halfint_pct = (100*(2.7478*spill_std)/(ite**0.5))/spill_mean
	effabs_halfint_pct = (100*(2.7478*effabs_std)/(ite**0.5))/effabs_mean
	e_per_ray = 37*218.0/heliohitsT
	H_norm = H*(e_per_ray/area)
	print("Halfbound% of spill, effabs and maxconc are "+str(spill_halfint_pct)+" "+str(effabs_halfint_pct)+" "+str(maxconc_halfint_pct))
	stop = timeit.default_timer()
	time_used = stop-start

	return [index,spill_mean,spill_halfint_pct,effabs_mean,effabs_halfint_pct,maxconc_mean,maxconc_halfint_pct,total_rays,time_used,H_norm,boundlist,extent]

def write_data(filename,lsofls):
	"""Writes data onto a given filename. lsofls refers to a list of list of results
	   such as [[0,1,2,3,4,5],[0,1,2,3,4,5]]. Normally the histogram and boundlist 
	   should be left out"""
	f = open(filename,'w')
	f.write("Index "+"Spill "+"Spill_bar "+"EffAbs "+"EffAbs_bar "+"MaxConc "+"MaxConc_bar "+"TotRays "+"TotTime "+"\n")
	for line in lsofls:	#take the first 6 values only
		f.write(str(line[0])+" "+str(line[1])+" "+str(line[2])+" "+str(line[3])+\
		" "+str(line[4])+" "+str(line[5])+" "+str(line[6])+" "+str(line[7])+" "+str(line[8])+"\n")
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
		
		
	
