#!/bin/python
# Copyright 2017 <Alessandro Erba>
import multiprocessing as mp
import os
import subprocess as sp
import numpy as np

root_dir = "."
solver = "../for-ch_solver"
num_cores = mp.cpu_count()
list_dataset = [f for f in os.listdir(root_dir) if f.endswith(".dat")]
checker = "pedibus_checker.py"
datasets = [300]

def start_solver_and_check_result(a, b, c, d, e, dataset, max_a, a_incr, max_b, b_incr, max_c, c_incr, max_d, d_incr, max_e, e_incr):	
	
	log = open("improvedParameters_"+ str(dataset) +".txt",'a+')

	if(a - max_a < 0):
		a_mod = a
	else: 
		a_mod = a - max_a
	if(b - max_b < 0):
		b_mod = b
	else:
		 b_mod = b - max_b
	if(c - max_c < 0):
		c_mod = 0
	else: 
		c_mod = c - max_c
	if(d - max_d < 0):
		d_mod = 0
	else: 
		d_mod = d - max_d
	if(e - max_e < 0):
		e_mod = 0
	else: 
		e_mod = e - max_e

	vet_a = np.linspace(a_mod, max_a, a_incr)
	vet_b = np.linspace(b_mod, max_b, b_incr)
	vet_c = np.linspace(c_mod, max_c, c_incr)
	vet_d = np.linspace(d_mod, max_d, d_incr)
	vet_e = np.linspace(e_mod, max_e, e_incr)
	
	for a_mod in vet_a:
		for b_mod in vet_b:
			for c_mod in vet_c:
				for d_mod in vet_d:
					for e_mod in vet_e:
						param = " -t 0 -a "+ str(a_mod) +" -b "+ str(b_mod) +" -c "+ str(c_mod) +" -d "+ str(d_mod) +" -e "+ str(e_mod)
						command = solver + param +" pedibus_" + str(dataset)  + ".dat " 
						print("Execution:", command)
						log.write("---------------------------------")
						log.write('\n' + param + '\n' + str(dataset) + '\n')
						log.flush()
						try:
						    import time
						    start_time = time.time()
						    sp.call(command, timeout = 900, shell = True)
						except sp.TimeoutExpired:
						    log.write("timeout expired")
						log.write("\n empleased time: %s \n" % (time.time() - start_time))    
						param_checker = " pedibus_" + str(dataset) + ".dat ../pedibus_" + str(dataset) + ".sol"						
						sp.call("python " + checker + param_checker, stdout = log, shell = True)
						log.flush()
						
						log.write("--------------------------------\n")


for dataset in datasets:
	if (dataset == 10):
		a = 0.1;
		b = 1;
		c = 100;
		d = 0.4;
		e = 0.1;
		#current_min_leaf = 444
		max_a = 0.2
		a_incr = 5
		max_b = 3
		b_incr = 7
		max_c = 200
		c_incr = 201
		max_d = 0.9
		d_incr = 101
		max_e = 0.9
		e_incr = 101
		start_solver_and_check_result(a, b, c, d, e, dataset, max_a, a_incr, max_b, b_incr, max_c, c_incr, max_d, d_incr, max_e, e_incr)

	if (dataset == 20):
		a = 0.1;
		b = 1;
		c = 100;
		d = 0.4;
		e = 0;
		max_a = 0.2
		a_incr = 4
		max_b = 3
		b_incr = 0.5
		max_c = 200
		c_incr = 200
		max_d = 0.9
		d_incr = 100
		max_e = 0.9
		e_incr = 100
		start_solver_and_check_result(a, b, c, d, e, dataset, max_a, a_incr, max_b, b_incr, max_c, c_incr, max_d, d_incr, max_e, e_incr)

	if (dataset == 30):
		a = 0.1;
		b = 1;
		c = 100;
		d = 0.4;
		e = 0.1;
		max_a = 0.2
		a_incr = 4
		max_b = 3
		b_incr = 0.5
		max_c = 200
		c_incr = 200
		max_d = 0.9
		d_incr = 100
		max_e = 0.9
		e_incr = 100
		start_solver_and_check_result(a, b, c, d, e, dataset, max_a, a_incr, max_b, b_incr, max_c, c_incr, max_d, d_incr, max_e, e_incr)

	if (dataset == 50):
		a = 0.1;
		b = 1;
		c = 100;
		d = 0.4;
		e = 0;
		max_a = 0.2
		a_incr = 5
		max_b = 3
		b_incr = 7
		max_c = 200
		c_incr = 21
		max_d = 0.9
		d_incr = 11
		max_e = 0.9
		e_incr = 11
		start_solver_and_check_result(a, b, c, d, e, dataset, max_a, a_incr, max_b, b_incr, max_c, c_incr, max_d, d_incr, max_e, e_incr)

	if (dataset == 80):
		a = 0.1;
		b = 1;
		c = 100;
		d = 0.4;
		e = 0.1;
		max_a = 0.2
		a_incr = 5
		max_b = 3
		b_incr = 7
		max_c = 200
		c_incr = 21
		max_d = 0.9
		d_incr = 11
		max_e = 0.9
		e_incr = 11
		start_solver_and_check_result(a, b, c, d, e, dataset, max_a, a_incr, max_b, b_incr, max_c, c_incr, max_d, d_incr, max_e, e_incr)


	if (dataset == 100):
		a = 0.1;
		b = 1;
		c = 100;
		d = 0.4;
		e = 0.1;
		max_a = 0.2
		a_incr = 5
		max_b = 3
		b_incr = 7
		max_c = 200
		c_incr = 21
		max_d = 0.9
		d_incr = 11
		max_e = 0.9
		e_incr = 11
		start_solver_and_check_result(a, b, c, d, e, dataset, max_a, a_incr, max_b, b_incr, max_c, c_incr, max_d, d_incr, max_e, e_incr)


	if (dataset == 150):
		a = 0;
		b = 1;
		c = 30;
		d = 0.4;
		e = 0.01;
		max_a = 0.2
		a_incr = 5
		max_b = 3
		b_incr = 7
		max_c = 200
		c_incr = 21
		max_d = 0.9
		d_incr = 11
		max_e = 0.9
		e_incr = 101
		start_solver_and_check_result(a, b, c, d, e, dataset, max_a, a_incr, max_b, b_incr, max_c, c_incr, max_d, d_incr, max_e, e_incr)


	if (dataset == 200):
		a = 0;
		b = 1;
		c = 30;
		d = 0.7;
		e = 0.01;
		max_a = 0.2
		a_incr = 5
		max_b = 3
		b_incr = 7
		max_c = 200
		c_incr = 21
		max_d = 0.9
		d_incr = 11
		max_e = 0.9
		e_incr = 101
		start_solver_and_check_result(a, b, c, d, e, dataset, max_a, a_incr, max_b, b_incr, max_c, c_incr, max_d, d_incr, max_e, e_incr)


	if (dataset == 250):
		a = 0.1;
		b = 1;
		c = 100;
		d = 0.7;
		e = 0;
		max_a = 0.2
		a_incr = 5
		max_b = 3
		b_incr = 7
		max_c = 200
		c_incr = 21
		max_d = 0.9
		d_incr = 11
		max_e = 0.9
		e_incr = 101
		start_solver_and_check_result(a, b, c, d, e, dataset, max_a, a_incr, max_b, b_incr, max_c, c_incr, max_d, d_incr, max_e, e_incr)

	if (dataset == 300):
		a = 0.1;
		b = 1;
		c = 100;
		d = 0.7;
		e = 0;
		max_a = 0.2
		a_incr = 5
		max_b = 3
		b_incr = 7
		max_c = 150
		c_incr = 11
		max_d = 0.9
		d_incr = 10
		max_e = 0.9
		e_incr = 100
		start_solver_and_check_result(a, b, c, d, e, dataset, max_a, a_incr, max_b, b_incr, max_c, c_incr, max_d, d_incr, max_e, e_incr)


