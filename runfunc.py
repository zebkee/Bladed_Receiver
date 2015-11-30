# Import necessary functions
from functions.zebfn import *
from models.receiver import*
import timeit

# Optimization functions
def bladed(A,B,W1,W2,H,D,R):
	"""Just generates a tower scene, returns the scene"""
	a,b,c = true_bladed_receiver(A,B,W1,W2,H,D,R)
	
	scene = TowerScene(reclist = a, recobj = b, objlist = c, z_offset = 6.1)
	scene.aim_field()
	return scene

def bladed_trace(index,A=0.9,B=4,W1=1,W2=1,H=1,D=0.9,R=1,tot_rays=1000000,ray_bund=100000):
	"""Does the tracing for tot_rays of rays using ray_bund rays at a time. Index is the value of the manipulated variable. Returns a list of [index,spill,effabs,impratio,max_conc]"""
	start = timeit.default_timer()
	n = 0
	energyT = 0
	subhitsT = 0

	scene = bladed(A,B,W1,W2,H,D,R)
	scene.trace(hstat_rays=int(ray_bund/218),render=False)
	H1, boundlist, extent, area = scene.z_histdata(50)
	n += 218*int(ray_bund/218)
	energy,subhits = scene.z_totalenergy(A)
	energyT += energy
	subhitsT += subhits
	binarea = area
	H2 = H1
	#print("Processed "+str(n)+" rays.")
	while n < tot_rays:
		scene = bladed(A,B,W1,W2,H,D,R)
		scene.trace(hstat_rays=int(ray_bund/218),render=False)
		H1, boundlist, extent, area = scene.z_histdata(50)
		n += 218*int(ray_bund/218)
		energy,subhits = scene.z_totalenergy(A)
		energyT += energy
		subhitsT += subhits
		H2 += H1
	print("Processed "+str(n)+" rays.")
	spill = 1.0*(n-subhitsT)/(n)
	effabs = energyT/(subhitsT)
	#effabs = energyT/((1-spill)*n)
	impratio = effabs/A
	e_per_ray = 218.0/n #this is in kWatts
	max_raw = N.amax(H2)
	factor = e_per_ray/binarea
	max_conc = max_raw*factor
	scene = None
	stop = timeit.default_timer()
	time = stop-start
	#print("MaxRaw = "+str(max_raw))
	

	return [index,spill,effabs,impratio,max_conc,n,time]

def bladed_traceMP(proc,index,A=0.9,B=4,W1=1,W2=1,H=1,D=0.9,R=1,tot_rays=1000000,ray_bund=100000):
	"""Does the tracing for tot_rays of rays using ray_bund rays at a time. Index is the value of the manipulated variable. Returns a list of [index,spill,effabs,impratio,max_conc]"""
	start = timeit.default_timer()
	n = 0
	energyT = 0
	subhitsT = 0
	#rppph is rays per heliostat per processor
	# New: rph is rays per heliostat
	rph = int(ray_bund/218)
	#print(rppph)

	scene = bladed(A,B,W1,W2,H,D,R)
	scene.traceMP(procs = proc,hstat_rays=rph,render=False)
	H1, boundlist, extent, area = scene.z_histdata(50)
	n += rph*proc*218
	energy,subhits = scene.z_totalenergy(A)
	energyT += energy
	subhitsT += subhits
	binarea = area
	H2 = H1
	scene = None
	#print("Processed "+str(n)+" rays.")
	while n < tot_rays:
		scene = bladed(A,B,W1,W2,H,D,R)
		scene.traceMP(procs=proc,hstat_rays=rph,render=False)
		H1, boundlist, extent, area = scene.z_histdata(50)
		n += rph*proc*218
		energy,subhits = scene.z_totalenergy(A)
		energyT += energy
		subhitsT += subhits
		H2 += H1
		scene = None
	print("Processed "+str(n)+" rays.")
	spill = 1.0*(n-subhitsT)/(n)
	effabs = energyT/(subhitsT)
	#effabs = energyT/((1-spill)*n)
	impratio = effabs/A
	e_per_ray = 218.0/n #this is in kWatts
	max_raw = N.amax(H2)
	factor = e_per_ray/binarea
	max_conc = max_raw*factor
	scene = None
	stop = timeit.default_timer()
	time = stop-start
	#print("MaxRaw = "+str(max_raw))
	

	return [index,spill,effabs,impratio,max_conc,n,time]

def write_data(filename,result):
	"""Writes theresults to a csv file defined by the input. result is an lsofls [[1,2,3,4,5],[1,2,3,4,5]]"""
	f = open(filename,'w')
	f.write("Index "+"Spill "+"EffAbs "+"ImpRatio "+"MaxConc "+"TotRays "+"TotTime "+"\n")
	for line in result:
		f.write(str(line[0])+" "+str(line[1])+" "+str(line[2])+" "+str(line[3])+" "+str(line[4])+" "+str(line[5])+" "+str(line[6])+"\n")
	f.close()

def bladed_count(A,B,W1,W2,H,D,R,nmax,delta):
	"""Plots the variation of results against the number of rays used to test for stability"""
	#Initialize the results table
	results = open("/home/zebedee/Desktop/result_bladed_count.csv",'w')
	results.write("Rays "+"Spill "+"EffAbs "+"Impratio "+"MaxConc"+"\n")

	n = 0
	energyT = 0
	subhitsT= 0

	scene = bladed(A,B,W1,W2,H,D,R)
	scene.trace(hstat_rays=int(delta/218),render=False)
	H1, boundlist, extent, area = scene.z_histdata(50)
	n += 218*int(delta/218)
	energy,subhits = scene.z_totalenergy(A)
	#print(energy)
	#print(subhits)
	energyT += energy
	subhitsT += subhits
	binarea = area
	spill = 1.0*(n-subhitsT)/(n)
	effabs = energyT/(subhitsT)
	impratio = effabs/A
	e_per_ray = 218.0/n #this is in kWatts
	factor = e_per_ray/binarea
	H2 = H1
	max_raw = N.amax(H2)
	#H_norm = H1*(factor)
	max_conc = max_raw*factor
	results.write(str(n)+" "+str(spill)+" "+str(effabs)+" "+str(impratio)+" "+str(max_conc)+"\n")
	#print("Processed "+str(n)+" rays.")

	while n < nmax:
		scene = bladed(A,B,W1,W2,H,D,R)
		scene.trace(hstat_rays=int(delta/218),render=False)
		H1, boundlist, extent, area = scene.z_histdata(50)
		n += 218*int(delta/218)
		energy,subhits = scene.z_totalenergy(A)
		#print(energy)
		#print(subhits)
		energyT += energy
		subhitsT += subhits
		spill = 1.0*(n-subhitsT)/(n)
		effabs = energyT/(subhitsT)
		impratio = effabs/A
		e_per_ray = 218.0/n #this is in kWatts
		factor = e_per_ray/binarea
		H2 += H1
		max_raw = N.amax(H2)
		#H_norm = H1*(factor)
		max_conc = max_raw*factor
		#print(max_conc)
	
		#Write the results file
		results.write(str(n)+" "+str(spill)+" "+str(effabs)+" "+str(impratio)+" "+str(max_conc)+"\n")
		print("Processed "+str(n)+" rays.")
	results.close()

def bladed_number(a=0.9,x=0,y=0,z=0,rx=0,ry=0,rz=0,w1=1,w2=1,h=1,d=0.9,r=1,bmin=3,bmax=10,filename="/home/zebedee/Desktop/result_bladed_number.csv"):
	"""Simulates bladed receivers with different number of blades using 10^6 rays each"""
	tot_rays = 1000000
	ray_bund = 100000
	lsofls = []
	B = bmin
	while B <= bmax:
		print("Blade number "+str(B))
		result = bladed_trace(B,a,B,w1,w2,h,d,r,tot_rays,ray_bund)
		lsofls.append(result)
		B += 1
	write_data(filename,lsofls)

def bladed_ratio(rmin=0.1,rmax=2.0,rint=0.1,a=0.9,x=0,y=0,z=0,rx=0,ry=0,w1=1,w2=1,d=0.9,b=8,filename="/home/zebedee/Desktop/result_bladed_ratio.csv"):
	"""Simulates bladed receivers with different thickness to spacing ratio using 10^6 rays each"""
	lsofls = []
	r = rmin
	while R <= rmax:
		print("Ratio of "+str(R))
		result = bladed_trace(r,a,B,w1,w2,h,d,r,1000000,100000)
		lsofls.append(result)
		R += rint
	write_data(filename,lsofls)
	
def compare(N,T,B):
	"""Compares the speed and results of linear and MP, does each one N number of times"""
	lsofls = []
	timelist = []
	n = 1
	c = 1
	#linear first
	#print("Linear "+str(c))
	while n <= 0:
		start = timeit.default_timer()
		result = bladed_trace(c,tot_rays=T,ray_bund=B)
		stop = timeit.default_timer()
		lsofls.append(result)
		timelist.append(stop-start)
		n += 1
	n = 1
	while c <= 3:
		print("MP "+str(c))
		while n <= N:
			start = timeit.default_timer()
			result = bladed_traceMP(c,c,tot_rays=T,ray_bund=B)
			stop = timeit.default_timer()
			lsofls.append(result)
			timelist.append(stop-start)
			n += 1
		n = 1
		c += 1
	print(timelist)

	write_data("/home/zebedee/Desktop/result_bladed_mp.csv",lsofls)

#Run code here
#start = timeit.default_timer()
#bladed_count(0.9,4,1.0,1.0,1.0,0.7,0.5,5000000,10000)
#bladed_number(bmin = 16, bmax = 17)
#stop = timeit.default_timer()

#runtime = (stop-start)/60 #in Minutes
#print("Runtime = "+str(runtime)+" minutes.")

#Compare Linear and MP time
compare(2,1000000,50000)
#start = timeit.default_timer()
#result = bladed_traceMP(3,3,tot_rays=1000000,ray_bund=50000)
#stop = timeit.default_timer()
#print(result)
#print(stop-start)
	
	
	
