from functions.blade_trace import *
from functions.tower_scene_mp_mark2 import *

def test_receiver(rays=2180,bins=50):
	result_ls = bladed_trace(1,rays,rays,bins,render_bool=True)
	n = 0
	while n < 6:
		print(result_ls[n])
		n += 1
	H_norm = result_ls[6]
	boundls = result_ls[7]
	extent = result_ls[8]
	show_hist(H_norm,boundls,1,extent)

def test_generic(rec_obj,surf_ls,crit_ls,absorp):
	""" Just test a generic receiver"""

	result_ls = generic_trace(1,rec_obj,surf_ls,crit_ls,absorp,rays_per_bund=218,tot_rays=100000,bins=50)
	n = 0
	while n < 6:
		print(result_ls[n])
		n += 1
	H_norm = result_ls[6]
	boundls = result_ls[7]
	extent = result_ls[8]
	show_hist(H_norm,boundls,1,extent)

def test_MP(rec_obj,surf_ls,crit_ls,absorp,procs=1):
    	rec_obj_copy = copy.deepcopy(rec_obj)
    	surf_ls_copy = copy.deepcopy(surf_ls)
   	crit_ls_copy = copy.deepcopy(crit_ls)
    	heliostat = "/home/zebedee/Tracer/Tracer/examples/sandia_hstat_coordinates.csv"
    	scene = TowerSceneMP(rec_obj_copy,surf_ls_copy,crit_ls_copy,heliostat,index=1,absorp=absorp)
    	result_ls = scene.traceMP(50,1000000,procs=procs,render=False,no_of_bins=50)

   	n = 0
   	while n < 6:
       		print(result_ls[n])
        	n += 1
    	H_norm = result_ls[6]
    	boundls = result_ls[7]
    	extent = result_ls[8]
    	show_hist(H_norm,boundls,1,extent)
    	print("Done")




width = 1
height = 1
rec_obj, surf_ls, crit_ls = cavity_aligned(width, height, absorptivity=0.8)
#test_generic(rec_obj, surf_ls, crit_ls,0.8)
test_MP(rec_obj,surf_ls,crit_ls,0.8,procs=3)
