from functions.tower_scene import *
from models.receiver import *
import matplotlib.pyplot as plt
import timeit

def generic_trace(index,rec_obj,surf_ls,crit_ls,absorp,rays_per_bund=10000,tot_rays=100000,bins=50,render_bool=True):
	"""Raytraces a generic receiver, where rec_obj is an assembled object, surf_ls is a list of all surfaces and
	   crit_ls is a list of all surfaces that are to be plotted in the histogram"""
	heliocsv = "/home/zebedee/Tracer/Tracer/examples/sandia_hstat_coordinates.csv"
	start = timeit.default_timer()
	
	A = absorp

	rays_used = 0
	hitsT = 0
	energyT = 0

	rph = int(rays_per_bund/218)

	scene = TowerScene(rec_obj=rec_obj,surf_ls = surf_ls, crit_ls = crit_ls, heliostat = heliocsv)
	scene.gen_plant()
	scene.aim_field()
	
	scene.trace(rph=rph,render = render_bool)
	H1, boundlist, extent, area = scene.hist_comb(no_of_bins=bins)
	energy,hits = scene.energies(A)
	energyT += energy
	hitsT += hits
	binarea = area #The first assignment is enough
	H2 = H1 #The first assignment is enough
	extent1 = extent
	boundlist1 = boundlist # First assignment is enough
	rays_used += rph*218
	del scene #Try to free up memory

	while rays_used < tot_rays:
		scene = TowerScene(rec_obj=rec_obj,surf_ls = surf_ls, crit_ls = crit_ls, heliostat = heliocsv)
		scene.gen_plant()
		scene.aim_field()
		scene.trace(rph=rph,render=False)
		H1, boundlist, extent, area = scene.hist_comb(no_of_bins=bins)
		energy, hits = scene.energies(A)
		energyT += energy
		hitsT += hits
		H2 += H1
		del scene #Try to free up memory
		rays_used += rph*218
	#Calculate the final results
	spill = 1.0*(rays_used-hitsT)/rays_used
	effabs = energyT/hitsT
	e_per_ray = 218.0/rays_used
	H_norm = H2*(e_per_ray/binarea)
	#print(H_norm)
	max_conc = N.amax(H_norm)
	stop = timeit.default_timer()
	time_used = stop-start

	return [index,spill,effabs,max_conc,rays_used,time_used,H_norm,boundlist1,extent1]
	
	

def bladed(absorp,no_of_blades,width_x,width_y,height,depth,ratio):
	"""Generates a scene using a bladed receiver and Sandia's heliostat"""
	# absorp = absorptivity of each receiver surface
	# width_x = width of the receiver in the x direction
	# width_y = width of the receiver in the y direction
	# height = height of the receiver
	# depth =  how deep the blade cavities are into the receiver
	# ratio = ratio of blade thickness to the spacing
	heliocsv = "/home/zebedee/Tracer/Tracer/examples/sandia_hstat_coordinates.csv"
	a,b,c = true_bladed_receiver(absorp,no_of_blades,width_x,width_y,height,depth,ratio)
	scene = TowerScene(rec_obj = a, surf_ls = b, crit_ls = c, heliostat=heliocsv)
	scene.gen_plant()
	scene.aim_field()
	return scene

def bladed_trace(index,rays_per_bund,tot_rays,bins,A=0.9,B=4,W1=1,W2=1,H=1,D=0.9,R=1,render_bool=False):
	"""Does raytracing for (tot_rays) rays using (rays_per_bund) rays at a time.
	   (index) is for the manipulated variable. Returns a list of results 
	   [index,spill,effabs,max_conc,rays_used,time_used,H_norm]"""
	# To access the H_norm, use results[6]
	start = timeit.default_timer()
	rays_used = 0
	hitsT = 0
	energyT = 0

	#Calculate the number of rays per heliostat
	rph = int(rays_per_bund/218)

	scene1 = bladed(A,B,W1,W2,H,D,R)
	scene1.trace(rph=rph,render = render_bool)
	H1, boundlist, extent, area = scene1.hist_comb(no_of_bins=bins)
	energy,hits = scene1.energies(A)
	energyT += energy
	hitsT += hits
	binarea = area #The first assignment is enough
	H2 = H1 #The first assignment is enough
	extent1 = extent
	boundlist1 = boundlist # First assignment is enough
	rays_used += rph*218
	del scene1 #Try to free up memory
	while rays_used < tot_rays:
		sceneN = bladed(A,B,W1,W2,H,D,R)
		sceneN.trace(rph=rph,render=False)
		H1, boundlist, extent, area = sceneN.hist_comb(no_of_bins=bins)
		energy, hits = sceneN.energies(A)
		energyT += energy
		hitsT += hits
		H2 += H1
		del sceneN #Try to free up memory
		rays_used += rph*218
	#Calculate the final results
	spill = 1.0*(rays_used-hitsT)/rays_used
	effabs = energyT/hitsT
	e_per_ray = 218.0/rays_used
	H_norm = H2*(e_per_ray/binarea)
	print(H_norm)
	max_conc = N.amax(H_norm)
	stop = timeit.default_timer()
	time_used = stop-start

	return [index,spill,effabs,max_conc,rays_used,time_used,H_norm,boundlist1,extent1]

def bladed_traceMP(index,rays_per_bund,tot_rays,bins,A=0.9,B=4,W1=1,W2=1,H=1,D=0.9,R=1,proc=1,render_bool=False):
	"""Does raytracing for (tot_rays) rays using (rays_per_bund) rays at a time.
	   (index) is for the manipulated variable. Returns a list of results 
	   [index,spill,effabs,max_conc,rays_used,time_used,H_norm,boundlist]"""
	
	start = timeit.default_timer()
	rays_used = 0
	hitsT = 0
	energyT = 0

	#Calculate the number of rays per heliostat
	rph = int(rays_per_bund/(218*proc))

	scene1 = bladed(A,B,W1,W2,H,D,R)
	scene1.traceMP(rph=rph,render = render_bool)
	H1, boundlist, extent, area = scene1.hist_comb(no_of_bins=bins)
	energy,hits = scene1.energies(A)
	energyT += energy
	hitsT += hits
	binarea = area #The first assignment is enough
	H2 = H1 #The first assignment is enough
	boundlist1 = boundlist # First assignment is enough
	extent1 = extent
	rays_used += rph*218*proc
	del scene1 #Try to free up memory
	while rays_used < tot_rays:
		sceneN = bladed(A,B,W1,W2,H,D,R)
		sceneN.traceMP(proc,rph=rph,render=False)
		H1, boundlist, extent, area = sceneN.hist_comb(no_of_bins=bins)
		energy, hits = sceneN.energies(A)
		energyT += energy
		hitsT += hits
		H2 += H1
		del sceneN #Try to free up memory
		rays_used += rph*218*proc
	#Calculate the final results
	spill = 1.0*(rays_used-hitsT)/rays_used
	effabs = energyT/hitsT
	e_per_ray = 218.0/rays_used
	H_norm = H2*(e_per_ray/binarea)
	max_conc = N.amax(H_norm)
	stop = timeit.default_timer()
	time_used = stop-start

	return [index,spill,effabs,max_conc,rays_used,time_used,H_norm,boundlist1,extent1]
	
def write_data(filename,lsofls):
	"""Writes data onto a given filename. lsofls refers to a list of list of results
	   such as [[0,1,2,3,4,5],[0,1,2,3,4,5]]. Normally the histogram and boundlist 
	   should be left out"""
	f = open(filename,'w')
	f.write("Index "+"Spill "+"EffAbs "+"MaxConc "+"TotRays "+"TotTime "+"\n")
	for line in result:	#take the first 6 values only
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
		
		
	
