'''
A tower/heliostat field example based on the Tower_gui.py example.
It has been cleaned from GUI and mayaVi stuff and given the CoIn rendering import and capacity.
'''
#import traits.api as t_api
#import traitsui.api as tui

from visual import *

from models.receiver import *
from functions.zebfn import *
from tracer.CoIn_rendering.rendering import *

import numpy as N
from scipy.constants import degree

from tracer.ray_bundle import RayBundle
from tracer.sources import pillbox_sunshape_directions
from tracer.assembly import Assembly
from tracer.spatial_geometry import roty, rotx, rotz, rotation_to_z
from tracer.tracer_engine import TracerEngine
from tracer.optics_callables import *

from tracer.models.one_sided_mirror import one_sided_receiver
from tracer.models.heliostat_field import HeliostatField, radial_stagger, solar_vector

# For the flux map:
import matplotlib.pyplot as plt

class TowerScene():
    # Location of the sun:
    sun_az = 80. # degrees from positive X-axis
    sun_elev = 45. # degrees from XY-plane

    
    def __init__(self):
        self.gen_plant() 
   
    def gen_rays(self):
        sun_vec = solar_vector(self.sun_az*degree, self.sun_elev*degree)
        rpos = (self.pos + sun_vec).T
        direct = N.tile(-sun_vec, (self.pos.shape[0], 1)).T
        rays = RayBundle(rpos, direct, energy=N.ones(self.pos.shape[0]))
        
        return rays
    
    def gen_plant(self):
        # import custom coordinate file
        self.pos = N.loadtxt("/home/zebedee/Tracer/Tracer/examples/sandia_hstat_coordinates.csv", delimiter=',')
        self.pos *= 0.1
        # set heliostat field characteristics: 0.52m*0.52m, abs = 0, aim_h = 61
        self.field = HeliostatField(self.pos, 6.09e-1, 6.09e-1, 0, 6.1,1e-3)

	#add offset properties
	self.x_offset = x_offset
	self.y_offset = y_offset
	self.z_offset = z_offset

        #self.reclist, recobj = one_sided_receiver(1.0, 1.0, 0.2)
        #self.reclist, recobj = cavity_receiver(1.0, 1.0, absorp)
	self.reclist, recobj = thin_bladed_receiver(1.0, 1.0, absorp, blades)
        #rec_trans = rotx(N.pi/-2) # originally N.pi/2, changed to minus rotx(N.pi/-2)
	#print(recobj)
	rec_trans = rotx(0.)
        rec_trans[2,3] = z_offset # height of the tower original 6.1
	rec_trans[1,3] = y_offset #y-offset
	rec_trans[0,3] = x_offset #x-offset
        recobj.set_transform(rec_trans)
	if vertical == True:
		rec_trans2 = rotz(N.pi/2)*roty(N.pi/2)
		rec_trans2[2,3] = 6.1
		recobj.set_transform(rec_trans2)
	
        self.plant = Assembly(objects = [recobj], subassemblies=[self.field])
	global receiver
	receiver = recobj
    
    def aim_field(self):
        self.field.aim_to_sun(self.sun_az*degree, self.sun_elev*degree)
    
    def trace(self):
        """Generate a flux map using much more rays than drawn"""
        # Generate a large ray bundle using a radial stagger much denser
        # than the field.
        sun_vec = solar_vector(self.sun_az*degree, self.sun_elev*degree)
        
        hstat_rays = 7
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
        e.ray_tracer(rays, 100, 1e-9, tree=True)
        e.minener = 1e-9 # default 1e-5

	# Render:
        #trace_scene = Renderer(e)
        #trace_scene.show_rays()

        # Initialise a histogram of hits:
        #energy, pts = self.rec1.get_optics_manager().get_all_hits()
        #x, y = self.rec1.global_to_local(pts)[:2]
        #rngx = 0.55 #0.5
        #rngy = 0.55 #0.5
        
        #bins = 55 #50
        #H, xbins, ybins = N.histogram2d(x, y, bins, \
            #range=([-rngx,rngx], [-rngy,rngy]), weights=energy)
        #print(H, xbins, ybins)

        #extent = [ybins[0], ybins[-1], xbins[-1], xbins[0]]
        #plt.imshow(H, extent=extent, interpolation='nearest')
        #plt.colorbar()
	#plt.title("front")
        #plt.show()
	count = 0
	totalenergy = 0
	for face in self.reclist:
		energy, pts = face.get_optics_manager().get_all_hits()
		print(energy)
		print(pts)
		#print(energy.sum())
		subtotalenergy = energy.sum()
		totalenergy += subtotalenergy
		x, y = face.global_to_local(pts)[:2]
		rngx = 0.50 #0.5
		rngy = 0.50 #0.5

		bins = 50 #50
		H, xbins, ybins = N.histogram2d(x,y,bins,range=([-rngx,rngx],[-rngy,rngy]), weights=energy)
		#print(H, xbins, ybins)
		extent = [ybins[0], ybins[-1], xbins[-1], xbins[0]]
        	plt.imshow(H, extent=extent, interpolation='nearest')
        	plt.colorbar()
		if count == 0:
			plt.title("Front")
		else:
			plt.title(str(count))
        	plt.show()
		count += 1
	print('total energy below')
	print(totalenergy)
	#assume 218 heliostats
	eff = totalenergy/(218.0*hstat_rays)
	print("Blade apsorptivity = "+str(absorp))
	print("Effective absorptivity = "+str(eff))

	#start the 3D vpythonhere
	#create a dxdxd dictionary of small cubes
	w = span
	d = divs
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
				#print(tuppy)
				cubes[tuply] = 0.0
				#print([x,y,z])
				z += 1
			y += 1
			z = 0
		x += 1
		y = 0
	print(len(cubes))

	# now to assign each hit to a small cube
	delta = w/(2.0*d - 2.0)
	print(delta)
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
				inside = abs(key[0] - x) <= delta and abs(key[1] - y) <= delta and abs(key[2] - z) <= delta
				#print(inside)
				if inside == True:
					cubes[key] += energy[n]
					#print energy[n]
			n += 1
	#for item in cubes.values():
		#print item
	print(sum(cubes.values()))

	#delete empty boxes
	for item in cubes.keys():
		if cubes[item] == 0.:
			del cubes[item]

	#now time to render it
	l = w/(d-1.0)
	for item in cubes.keys():
		if cubes[item] >= E:
			box(pos=vector(item), color=(1.0,0.0,0.1), size=(l,l,l))
		elif cubes[item] > 0.75*E:
			box(pos=vector(item), color=(0.8,0.0,0.1), size=(l,l,l))
		elif cubes[item] > 0.5*E:
			box(pos=vector(item), color=(0.6,0.0,0.1), size=(l,l,l))
		elif cubes[item] > 0.25*E:
			box(pos=vector(item), color=(0.4,0.0,0.1), size=(l,l,l))
		elif cubes[item] > 0.0*E:
			box(pos=vector(item), color=(0.2,0.0,0.1), size=(l,l,l))	
	
#set parameters here
blades = 5
absorp = 0.7
vertical = False
x_offset = 0.
y_offset = 0.
z_offset = 6.1
span = 1.2
divs = 40.
E = 1.4

scene = TowerScene()
scene.aim_field()
scene.trace()
