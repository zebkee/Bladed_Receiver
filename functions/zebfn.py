#contains analytical functions edited by Zeb

# To use, refer to as from functions.zebfn import *

# Import Python Packages
from visual import *
import numpy as N
from numpy import r_
from scipy.constants import degree
import matplotlib.pyplot as plt
# Import Tracer Packages
from tracer.ray_bundle import RayBundle
from tracer.sources import pillbox_sunshape_directions
from tracer.assembly import Assembly
from tracer.spatial_geometry import roty, rotx, rotz, rotation_to_z
from tracer.tracer_engine import TracerEngine
from tracer.optics_callables import *
# Import Tracer Models
from tracer.models.heliostat_field import HeliostatField, radial_stagger, solar_vector
# Import Zeb Models
from models.receiver import *
# Import Coin3D Renderer
from tracer.CoIn_rendering.rendering import *

def z_hist1d(H,binny=50):
   	"""Creates a 1d histogram to disply flux uniformity"""
	bigarray = [item for sublist in H for item in sublist]
	maxconc = max(bigarray)
	print("The maximum concentration is "+str(maxconc)+" W per square metre")
	H_1, binedge = N.histogram(bigarray,bins = binny)
	plt.hist(H_1,bins = binedge)
	plt.show()
	#print(N.std(bigarray))

# Creates a scene of the tower
class TowerScene():
    # reclist and recobj are the two objects returned from thin_bladed_rec(1.0,1.0,absorp,blades)
    # Location of the sun


    def __init__(self,reclist,recobj,objlist,sun_azi = 80.,sun_eleva = 45.,x_offset = 0.,y_offset = 0.,z_offset = 6.1):
	self.sun_az = sun_azi # degrees from positive X-axis
        self.sun_elev = sun_eleva # degrees from XY-plane
        self.reclist = reclist
        self.recobj = recobj
	self.objlist = objlist
        #add offset properties
        self.x_offset = x_offset
        self.y_offset = y_offset
        self.z_offset = z_offset
        self.gen_plant() 
   
    def gen_rays(self):
        sun_vec = solar_vector(self.sun_az*degree, self.sun_elev*degree)
        rpos = (self.pos + sun_vec).T
        direct = N.tile(-sun_vec, (self.pos.shape[0], 1)).T
        rays = RayBundle(rpos, direct, energy=N.ones(self.pos.shape[0]))
        
        return rays
    
    def gen_plant(self):
	# reclist and recobj are the two objects returned from thin_bladed_rec(1.0,1.0,absorp,blades)
        # import custom coordinate file from my folder
        self.pos = N.loadtxt("/home/zebedee/Tracer/Tracer/examples/sandia_hstat_coordinates.csv", delimiter=',')
        self.pos *= 0.1
        # set heliostat field characteristics: 0.52m*0.52m, abs = 0, aim_h = 61
        self.field = HeliostatField(self.pos, 6.09e-1, 6.09e-1, 0, 6.1,1e-3)

        #self.reclist, recobj = one_sided_receiver(1.0, 1.0, 0.2)
        #self.reclist, recobj = cavity_receiver(1.0, 1.0, absorp)
	#self.reclist, recobj = thin_bladed_receiver(1.0, 1.0, absorp, blades)
        #rec_trans = rotx(N.pi/-2) # originally N.pi/2, changed to minus rotx(N.pi/-2)
	#print(recobj)
	rec_trans = rotx(0.)
        rec_trans[2,3] = self.z_offset # height of the tower original 6.1
	rec_trans[1,3] = self.y_offset #y-offset
	rec_trans[0,3] = self.x_offset #x-offset
        self.recobj.set_transform(rec_trans)
	#if vertical == True:
		#rec_trans2 = rotz(N.pi/2)*roty(N.pi/2)
		#rec_trans2[2,3] = 6.1
		#self.recobj.set_transform(rec_trans2)
	
        self.plant = Assembly(objects = [self.recobj], subassemblies=[self.field])
    
    def aim_field(self):
        self.field.aim_to_sun(self.sun_az*degree, self.sun_elev*degree)
    
    def trace(self, hstat_rays = 5, iters = 100, minE = 1e-5, render=True):
        """Generate a flux map using much more rays than drawn"""
        # Generate a large ray bundle using a radial stagger much denser
        # than the field.
	# hstat_rays refers to number of rays hitting each heliostat mirror
	# iters is no. of iterations, originally 100
	# minE refers to minimum energy of a ray before it gets discarded org=1e-9
	self.hstat_rays = hstat_rays
        sun_vec = solar_vector(self.sun_az*degree, self.sun_elev*degree)
        
        num_rays = hstat_rays*len(self.field.get_heliostats())
        rot_sun = rotation_to_z(-sun_vec)
        direct = N.dot(rot_sun, pillbox_sunshape_directions(num_rays, 0.00465))
        
        xy = N.random.uniform(low=-0.25, high=0.25, size=(2, num_rays))
        base_pos = N.tile(self.pos, (hstat_rays, 1)).T
        base_pos += N.dot(rot_sun[:,:2], xy)
        
        base_pos -= direct
        rays = RayBundle(base_pos, direct, energy=N.ones(num_rays))
        
        # Perform the trace:
        e = TracerEngine(self.plant)
        e.ray_tracer(rays, iters, minE, tree=True)
        e.minener = minE # default 1e-5

	if render == True:
        	trace_scene = Renderer(e)
        	trace_scene.show_rays()

    def z_hist(self,no_of_bins=50):
	""" Generates a histogram of all surfaces present within self.reclist and displays them one after another when  manually closed
	    Also prints the total energy hits of all surfaces"""
	count = 0
	totalenergy = 0
	for face in self.reclist:
		energy, pts = face.get_optics_manager().get_all_hits()
		#print(face.mesh(1))
		#calculate dimensions of each plate
		zarray = face.mesh(1)
		print(zarray)
		x_1 = zarray[0][0][0]
		x_2 = zarray[0][0][1]
		x_3 = zarray[0][1][0]
		y_1 = zarray[1][0][0]
		y_2 = zarray[1][0][1]
		y_3 = zarray[1][1][0]
		z_1 = zarray[2][0][0]
		z_2 = zarray[2][0][1]
		z_3 = zarray[2][1][0]
		w = ((x_2-x_1)**2+(y_2-y_1)**2+(z_2-z_1)**2)**0.5
		h = ((x_3-x_1)**2+(y_3-y_1)**2+(z_3-z_1)**2)**0.5
		print(w)
		print(h)
		#print(face.rotation)
		#print(pts)
		#print(energy.sum())
		subtotalenergy = energy.sum()
		totalenergy += subtotalenergy
		y, x = face.global_to_local(pts)[:2]
		#print(x)
		#print(y)
		rngx = h/2 #0.5
		rngy = w/2 #0.5

		bins = no_of_bins #50
		H, xbins, ybins = N.histogram2d(x,y,bins,range=([-rngx,rngx],[-rngy,rngy]), weights=energy)
		#print(H, xbins, ybins)
		extent = [ybins[0], ybins[-1], xbins[0], xbins[-1]]
        	plt.imshow(H, extent=extent, interpolation='nearest')
        	plt.colorbar()
		if count == 0:
			plt.title("Front")
		else:
			plt.title(str(count))
        	plt.show()
		count += 1
	print('The total energy hit is :'+str(totalenergy))

    def z_histall(self,no_of_bins=100):
	""" Generates a histogram of all surfaces present within self.reclist and displays them one after another when  manually closed
	    Also prints the total energy hits of all surfaces"""
	count = 0
	totalenergy = 0
	X_offset = 0
	all_X = []
	all_Y = []
	all_E = []
	boundlist = [0]
	for face in self.reclist:
		#Make empty list for boundary lines
		energy, pts = face.get_optics_manager().get_all_hits()
		#print(face.mesh(1))
		#calculate dimensions of each plate
		#get corners of the plate
		zarray = face.mesh(1)
		#Bottom Left Corner	1
		#Bottom Right Corner	2
		#Top Left Corner	3
		x_1 = zarray[0][1][1]
		x_2 = zarray[0][0][1]
		x_3 = zarray[0][1][0]
		y_1 = zarray[1][1][1]
		y_2 = zarray[1][0][1]
		y_3 = zarray[1][1][0]
		z_1 = zarray[2][1][1]
		z_2 = zarray[2][0][1]
		z_3 = zarray[2][1][0]
		BLC = r_[x_1,y_1,z_1]
		BRC = r_[x_2,y_2,z_2]
		TLC = r_[x_3,y_3,z_3]
		#print(BLC)
		#print(BRC)
		#print(TLC)
		#calculate width and height
		w = ((x_2-x_1)**2+(y_2-y_1)**2+(z_2-z_1)**2)**0.5
		h = ((x_3-x_1)**2+(y_3-y_1)**2+(z_3-z_1)**2)**0.5
		#get unit vectors of each side
		u = (r_[x_2,y_2,z_2]-r_[x_1,y_1,z_1])*(1.0/w)
		v = (r_[x_3,y_3,z_3]-r_[x_1,y_1,z_1])*(1.0/h)
		#get local points using dot product by unit vector
		#shifts the x position by width of previous combined plates
		local_X = list(((pts[0] - x_1)* u[0] + (pts[1] -y_1)* u[1] + (pts[2] -z_1) * u[2])+X_offset)
		local_Y = list((pts[0] - x_1)* v[0] + (pts[1] -y_1)* v[1] + (pts[2] -z_1) * v[2])
		#add them to a large list
		all_X += local_X
		all_Y += local_Y
		all_E += list(energy)
		subtotalenergy = energy.sum()
		totalenergy += subtotalenergy
		#print(x)
		#print(y)
		X_offset += w
		boundlist.append(X_offset)
	#bins = [int(no_of_bins*X_offset),no_of_bins]
	#print(len(all_E))
	#print(len(all_X))
	#print(len(all_Y))
	rngy = h
	rngx = X_offset
	bins = [no_of_bins,int(no_of_bins*X_offset)]

	#bins = [int(no_of_bins*X_offset),no_of_bins] #50
	H, ybins, xbins = N.histogram2d(all_Y,all_X,bins,range=([0,rngy],[0,rngx]), weights=all_E)
	#H, xbins, ybins = N.histogram2d(N.asarray(all_Y),N.asarray(all_X),bins,range=([0,h],[0,X_offset]), weights=N.asarray(all_E))
	#print(H, xbins, ybins)
	#print(H)
       	extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]
        img = plt.imshow(H, extent=extent, interpolation='nearest')
        plt.colorbar()
	#plt.vlines(boundlist,0,h,colors='k')
	y = r_[0,h]
	n = 1
	for bound in boundlist:
		globals()['line%s' % n] = plt.plot(r_[bound,bound],y,color='k')
		#print(globals()['line%s' % n])
		n += 1
	plt.xlim(0,X_offset)
        plt.show()
	return H

    def z_totalenergy(self,absorb):
	"""Returns the total energy absorbed by the surface, basicall the sum of all hits"""
	totalenergy = 0
	hits = 0
	for face in self.objlist:
		energy, pts = face.get_optics_manager().get_all_hits()
		totalenergy += energy.sum()
		#print(list(energy))
		hits += len([i for i in list(energy) if i > 0.999*absorb])
	return totalenergy, hits

    def z_histdata(self,no_of_bins=100):
	"""Generates just the basic results of a histogram array and a list of plate boundaries"""
	count = 0
	totalenergy = 0
	X_offset = 0
	all_X = []
	all_Y = []
	all_E = []
	boundlist = [0]
	for face in self.reclist:
		#Make empty list for boundary lines
		energy, pts = face.get_optics_manager().get_all_hits()
		#print(face.mesh(1))
		#calculate dimensions of each plate
		#get corners of the plate
		zarray = face.mesh(1)
		#Bottom Left Corner	1
		#Bottom Right Corner	2
		#Top Left Corner	3
		x_1 = zarray[0][1][1]
		x_2 = zarray[0][0][1]
		x_3 = zarray[0][1][0]
		y_1 = zarray[1][1][1]
		y_2 = zarray[1][0][1]
		y_3 = zarray[1][1][0]
		z_1 = zarray[2][1][1]
		z_2 = zarray[2][0][1]
		z_3 = zarray[2][1][0]
		BLC = r_[x_1,y_1,z_1]
		BRC = r_[x_2,y_2,z_2]
		TLC = r_[x_3,y_3,z_3]
		#calculate width and height
		w = ((x_2-x_1)**2+(y_2-y_1)**2+(z_2-z_1)**2)**0.5
		h = ((x_3-x_1)**2+(y_3-y_1)**2+(z_3-z_1)**2)**0.5
		#get unit vectors of each side
		u = (r_[x_2,y_2,z_2]-r_[x_1,y_1,z_1])*(1.0/w)
		v = (r_[x_3,y_3,z_3]-r_[x_1,y_1,z_1])*(1.0/h)
		#get local points using dot product by unit vector
		#shifts the x position by width of previous combined plates
		local_X = list(((pts[0] - x_1)* u[0] + (pts[1] -y_1)* u[1] + (pts[2] -z_1) * u[2])+X_offset)
		local_Y = list((pts[0] - x_1)* v[0] + (pts[1] -y_1)* v[1] + (pts[2] -z_1) * v[2])
		#add them to a large list
		all_X += local_X
		all_Y += local_Y
		all_E += list(energy)
		subtotalenergy = energy.sum()
		totalenergy += subtotalenergy
		X_offset += w
		boundlist.append(X_offset)
	rngy = h
	rngx = X_offset
	bins = [no_of_bins,int(no_of_bins*X_offset)]
	H, ybins, xbins = N.histogram2d(all_Y,all_X,bins,range=([0,rngy],[0,rngx]), weights=all_E)
       	extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]
	binarea = (float(h)/no_of_bins)*(float(X_offset)/int(no_of_bins*X_offset))
	return H, boundlist, extent, binarea

    def z_box_approx(self,span,divs):
	""" Generates a dictionary of tuple_coordinates to energy_floats, and box_side_length according to the span of approximation, usuall equal to the max length of sides of receiver and the number of divisions, bigger number gives higher resolution """
	#generates the tuple coordinates with zero energy to start with
	w = span
	d = divs
	x_offset = self.x_offset
	y_offset = self.y_offset
	z_offset = self.z_offset
	subint = w/(d-1)
	cubes = {}
	x = 0
	y = 0
	z = 0
	while x < d:
		while y < d:
			while z < d: #sets up coordinates fo each small cube
				x_n = x*subint - 0.5*w + x_offset
				y_n = y*subint - 0.5*w + y_offset
				z_n = z*subint - 0.5*w + z_offset
				tuply = (x_n,y_n,z_n)
				cubes[tuply] = 0.0
				z += 1
			y += 1
			z = 0
		x += 1
		y = 0
	print("Approximation done using "+str(len(cubes))+" small boxes")

	#Time to assign each hit to its correct box
	delta = w/(2.0*d - 2.0)
	for surfs in self.reclist:
		energy, pts = surfs.get_optics_manager().get_all_hits()
		hitno = len(energy)
		#print(len(energy))
		n = 0
		while n < hitno:
			x = pts[0][n]
			y = pts[1][n]
			z = pts[2][n]
			for key in cubes.keys():
				#print(cubes[key])
				#print(key[0])
				if abs(key[0] - x) <= delta and abs(key[1] - y) <= delta and abs(key[2] - z) <= delta:
					cubes[key] += energy[n]
			n += 1
	print("Total energy from approximation method is: "+str(sum(cubes.values())))

	#Delete irrelevant boxes
	newcube = {}
	for item in cubes.keys():
		if cubes[item] > 0.:
			newcube[item] = cubes[item]
	newno = len(newcube)
	perc = 100.0*len(newcube)/len(cubes)
	print("A total of "+str(newno)+" boxes were used giving a percentage of "+ str(perc) +"%")
	l = w/(d-1.0)
	self.cube = newcube
	self.l = l

    def z_box_approx_2(self):
	""" Generates a dictionary of tuple_coordinates to energy_floats, and box_side_length according to the span of approximation, usuall equal to the max length of sides of receiver and the number of divisions, bigger number gives higher resolution """
	#generates the tuple coordinates with zero energy to start with
	cubes = {}
	for surfs in self.reclist:
		mesh = surfs.mesh(100)
		m = 1
		n = 1
		while n < 100:
			while m < 100:
				x = mesh[0][n][m]
				y = mesh[1][n][m]
				z = mesh[2][n][m]
				cubes[(x,y,z)] = 0.0
				m += 2
			m = 1
			n += 2
		#print(cubes)
		energy, pts = surfs.get_optics_manager().get_all_hits()
		hitno = len(energy)
		#print(len(energy))
		delta = 0.01
		n = 0
		while n < hitno:
			x = pts[0][n]
			y = pts[1][n]
			z = pts[2][n]
			for key in cubes.keys():
				#print(cubes[key])
				#print(key[0])
				if abs(key[0] - x) <= delta and abs(key[1] - y) <= delta and abs(key[2] - z) <= delta:
					cubes[key] += energy[n]
			n += 1
	print("Total energy from approximation method is: "+str(sum(cubes.values())))
	self.cube = cubes

    def z_renderV(self,E):
	"""Function that takes a dictionary and renders the boxes using Visual Python. E refers to maximum expected energy and l represents the side length of each box"""
	l = 0.02
	for item in self.cube.keys():
		if self.cube[item] >= E:
			box(pos=vector(item), color=(1.0,0.0,0.1), size=(l,l,l))
		elif self.cube[item] > 0.75*E:
			box(pos=vector(item), color=(0.8,0.0,0.1), size=(l,l,l))
		elif self.cube[item] > 0.5*E:
			box(pos=vector(item), color=(0.6,0.0,0.1), size=(l,l,l))
		elif self.cube[item] > 0.25*E:
			box(pos=vector(item), color=(0.4,0.0,0.1), size=(l,l,l))
		elif self.cube[item] > 0.0*E:
			box(pos=vector(item), color=(0.2,0.0,0.1), size=(l,l,l))
		else:
			box(pos=vector(item), color=(0.0,0.0,0.1), size=(l,l,l))
	scene.autoscale = False
	scene.userzoom = True
	scene.userspin = True






		
		

	
