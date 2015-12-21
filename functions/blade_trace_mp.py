from functions.tower_scene_mp import *
from functions.blade_trace import *
from models.receiver import *
import matplotlib.pyplot as plt
import timeit

import copy


def generic_trace_mp(index,rec_obj,surf_ls,crit_ls,absorp,rays_per_bund=10000,tot_rays=100000,bins=50,render_bool=True,T=[0.,0.,6.1,0.,0.,0.],sun_az = 0.,sun_elev = 34.9,procs=1):
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

	scene = TowerSceneMP(rec_obj=rec_obj_copy,surf_ls = surf_ls_copy, crit_ls = crit_ls_copy, heliostat = heliocsv,dx=T[0],dy=T[1],dz=T[2],rx=T[3],ry=T[4],rz=T[5],sun_az=sun_az,sun_elev=sun_elev)
	scene.aim_field()
	
	rpr = scene.traceMP(rays_per_bund,render = render_bool,procs=procs)
	H2, boundlist, extent, area = scene.hist_comb(no_of_bins=bins)
	print(N.sum(H2))
	energy,hits,heliohits = scene.energies(A)
	energyT = energy
	hitsT = hits
	heliohitsT = heliohits
	binarea = area #The first assignment is enough
	#H2 = H1 #The first assignment is enough
	extent1 = extent
	boundlist1 = boundlist # First assignment is enough
	rays_used += rpr
	del scene #Try to free up memory

	while rays_used < tot_rays:
		rec_obj_copy = copy.deepcopy(rec_obj)
        	surf_ls_copy = copy.deepcopy(surf_ls)
        	crit_ls_copy = copy.deepcopy(crit_ls)

		scene = TowerSceneMP(rec_obj=rec_obj_copy,surf_ls = surf_ls_copy, crit_ls = crit_ls_copy, heliostat = heliocsv)
		scene.aim_field()
		rpr = scene.traceMP(rays_per_bund,render=False,procs=procs)
		rays_used += rpr
		H1, boundlist, extent, area = scene.hist_comb(no_of_bins=bins)
		print(N.sum(H1))
		energy, hits, heliohits = scene.energies(A)
		energyT += energy
		hitsT += hits
		heliohitsT += heliohits
		H2 += H1
		del scene #Try to free up memory
	#Calculate the final results
	spill = 1.0*(rays_used-hitsT)/rays_used
	effabs = energyT/hitsT
	e_per_ray = 37*218.0/rays_used #each heliostat is 37m^2
	H_norm = H2*(e_per_ray/binarea)
	#print(H_norm)
	max_conc = N.amax(H_norm)
	stop = timeit.default_timer()
	time_used = stop-start

	return [index,spill,effabs,max_conc,rays_used,time_used,H_norm,boundlist1,extent1]

def precision_traceMP(index,rec_obj,surf_ls,crit_ls,rays_per_iter=200000,bins=100,precision=2.0,render_bool=False,T=[0.,0.,6.1,0.,0.,0.],sun_az = 0.,sun_elev = 34.9,procs=1):
	#precision is in %
	heliocsv = "/home/zebedee-u5277975/Tracer/Tracer/examples/sandia_hstat_coordinates.csv"
	
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
		scene = TowerSceneMP(rec_obj=rec_obj_copy,surf_ls = surf_ls_copy, crit_ls = crit_ls_copy,\
		        heliostat=heliocsv,dx=T[0],dy=T[1],dz=T[2],rx=T[3],ry=T[4],rz=T[5],sun_az=sun_az,sun_elev=sun_elev)
		scene.aim_field()
		rpr = scene.traceMP(rays_per_run, procs = procs, render = render_bool)
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
		print("Maxconc = "+str(maxconc_mean)+" +/- "+str(maxconc_std))
		maxconc_halfint_pct = (100*(2.7478*maxconc_std)/(ite**0.5))/maxconc_mean
		print("Half interval is now "+str(maxconc_halfint_pct)+"%")    
		spill_mean,spill_std = stats(spill_array)
		print("Spill = "+str(spill_mean)+" +/- "+str(spill_std))
		effabs_mean,effabs_std = stats(effabs_array)
		print("EffAbs = "+str(effabs_mean)+" +/- "+str(effabs_std))
		total_rays += rpr
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


