from functions.blade_trace import *

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

	result_ls = generic_trace(1,rec_obj,surf_ls,crit_ls,absorp,rays_per_bund=218,tot_rays=10000,bins=50,sun_az = 0.,sun_elev = 34.9)
	n = 0
	while n < 6:
		print(result_ls[n])
		n += 1
	H_norm = result_ls[6]
	boundls = result_ls[7]
	extent = result_ls[8]
	show_hist(H_norm,boundls,1,extent)

def test_genericMP(rec_obj,surf_ls,crit_ls,absorp,procs=1):
	""" Just test a generic receiver"""

	result_ls = single_traceMP(1,rec_obj,surf_ls,crit_ls,absorp,procs=procs,rays_per_bund=218,tot_rays=10000,bins=50)
	n = 0
	while n < 6:
		print(result_ls[n])
		n += 1
	H_norm = result_ls[6]
	boundls = result_ls[7]
	extent = result_ls[8]
	show_hist(H_norm,boundls,1,extent)




width = 1
height = 1
rec_obj, surf_ls, crit_ls = cavity_aligned(width, height, absorptivity=0.8)
test_generic(rec_obj, surf_ls, crit_ls,0.8)
