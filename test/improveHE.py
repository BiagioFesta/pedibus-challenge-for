#!/bin/python
# Copyright 2017 <Alessandro Erba>
import multiprocessing as mp
import os
import subprocess as sp
import numpy as np
import signal

num_cores = mp.cpu_count()
root_dir = "."
solver = "../for-ch_solver"
num_cores = mp.cpu_count()
list_dataset = [f for f in os.listdir(root_dir) if f.endswith(".dat")]
checker = "pedibus_checker.py"
datasets = [300, 250, 200, 150, 100, 80, 50]
print("Num cores detected: " + str(num_cores))
print("Dataset files detected:", list_dataset)
print("Solver detected: " + solver)
def start_solver_and_check_result(a, b, c, d, e, dataset, max_a, min_a, a_incr, max_b, min_b, b_incr, max_c, min_c, c_incr, max_d, min_d, d_incr, max_e, min_e, e_incr):	
	
	log = open("improvedParameters_"+ str(dataset) +".txt",'a+')

	vet_a = np.linspace(min_a, max_a, a_incr)
	vet_b = np.linspace(min_b, max_b, b_incr)
	vet_c = np.linspace(min_c, max_c, c_incr)
	vet_d = np.linspace(min_d, max_d, d_incr)
	vet_e = np.linspace(min_e, max_e, e_incr)
	print(vet_a)
	print(vet_b)
	print(vet_c)
	print(vet_d)
	print(vet_e)
	a_mod = min_a
	b_mod = min_b
	c_mod = min_c
	d_mod = min_d
	e_mod = min_e

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
						import time
						start_time = time.time()
						try:
							p = sp.Popen(command, shell = True, preexec_fn = os.setsid)
							output, error = p.communicate(timeout=900)
							log.write("\n empleased time: %s \n" % (time.time() - start_time))    
							param_checker = " pedibus_" + str(dataset) + ".dat ../pedibus_" + str(dataset) + ".sol"						
							process = sp.call("python " + checker + param_checker, stdout = log, shell = True)
							log.flush()
						except sp.TimeoutExpired:
							log.write("timeout expired")
							os.killpg(os.getpgid(p.pid), signal.SIGTERM)
						log.write("--------------------------------\n")

def launch_thread(dataset):
	
		if (dataset == 10):
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
			c_incr = 201
			max_d = 0.9
			d_incr = 101
			max_e = 0.9
			e_incr = 101
			min_a = 0
			min_b = 0
			min_c = 0
			min_d = 0
			min_e = 0
			start_solver_and_check_result(a, b, c, d, e, dataset, max_a, min_a, a_incr, max_b, min_b, b_incr, max_c, min_c, c_incr, max_d, min_d, d_incr, max_e, min_e, e_incr)

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
			min_a = 0
			min_b = 0
			min_c = 0
			min_d = 0
			min_e = 0
			start_solver_and_check_result(a, b, c, d, e, dataset, max_a, min_a, a_incr, max_b, min_b, b_incr, max_c, min_c, c_incr, max_d, min_d, d_incr, max_e, min_e, e_incr)

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
			min_a = 0
			min_b = 0
			min_c = 0
			min_d = 0
			min_e = 0
			start_solver_and_check_result(a, b, c, d, e, dataset, max_a, min_a, a_incr, max_b, min_b, b_incr, max_c, min_c, c_incr, max_d, min_d, d_incr, max_e, min_e, e_incr)

		if (dataset == 50):
			a = 0.1;
			b = 1;
			c = 100;
			d = 0.4;
			e = 0;
			max_a = 0.2
			min_a = 0
			a_incr = 3
			max_b = 3
			min_b = 0
			b_incr = 4
			max_c = 150
			min_c = 50
			c_incr = 11
			max_d = 0.9
			min_d = 0
			d_incr = 20
			max_e = 0.9
			min_e = 0
			e_incr = 20
			start_solver_and_check_result(a, b, c, d, e, dataset, max_a, min_a, a_incr, max_b, min_b, b_incr, max_c, min_c, c_incr, max_d, min_d, d_incr, max_e, min_e, e_incr)

		if (dataset == 80):
			a = 0.1;
			b = 1;
			c = 100;
			d = 0.4;
			e = 0.1;
			max_a = 0.2
			min_a = 0
			a_incr = 3
			max_b = 3
			min_b = 0
			b_incr = 4
			max_c = 150
			min_c = 50
			c_incr = 11
			max_d = 0.9
			min_d = 0
			d_incr = 20
			max_e = 0.9
			min_e = 0
			e_incr = 20
			start_solver_and_check_result(a, b, c, d, e, dataset, max_a, min_a, a_incr, max_b, min_b, b_incr, max_c, min_c, c_incr, max_d, min_d, d_incr, max_e, min_e, e_incr)


		if (dataset == 100):
			a = 0.1;
			b = 1;
			c = 100;
			d = 0.4;
			e = 0.1;
			max_a = 0.2
			min_a = 0
			a_incr = 3
			max_b = 3
			min_b = 0
			b_incr = 4
			max_c = 150
			min_c = 50
			c_incr = 11
			max_d = 0.9
			min_d = 0
			d_incr = 5.5
			max_e = 0.9
			min_e = 0
			e_incr = 5.5
			start_solver_and_check_result(a, b, c, d, e, dataset, max_a, min_a, a_incr, max_b, min_b, b_incr, max_c, min_c, c_incr, max_d, min_d, d_incr, max_e, min_e, e_incr)


		if (dataset == 150):
			a = 0;
			b = 1;
			c = 30;
			d = 0.4;
			e = 0.01;
			max_a = 0.3
			min_a = 0.1
			a_incr = 3
			max_b = 3
			min_b = 0
			b_incr = 4
			max_c = 60
			min_c = 20
			c_incr = 7.8
			max_d = 0.9
			min_d = 0
			d_incr = 5.5
			max_e = 0.9
			e_incr = 5.5
			start_solver_and_check_result(a, b, c, d, e, dataset, max_a, min_a, a_incr, max_b, min_b, b_incr, max_c, min_c, c_incr, max_d, min_d, d_incr, max_e, min_e, e_incr)


		if (dataset == 200):
			a = 0;
			b = 1;
			c = 30;
			d = 0.7;
			e = 0.01;
			max_a = 0.3
			min_a = 0.1
			a_incr = 3
			max_b = 3
			min_b = 1
			b_incr = 3	
			max_c = 60
			min_c = 20
			c_incr = 7.8
			max_d = 0.9
			min_d = 0
			d_incr = 5.5
			max_e = 0.9
			min_e = 0
			e_incr = 5.5
			start_solver_and_check_result(a, b, c, d, e, dataset, max_a, min_a, a_incr, max_b, min_b, b_incr, max_c, min_c, c_incr, max_d, min_d, d_incr, max_e, min_e, e_incr)


		if (dataset == 250): 
			a = 0.1;
			b = 1;
			c = 100;
			d = 0.7;
			e = 0;
			max_a = 0.3
			min_a = 0.1
			a_incr = 3
			max_b = 3
			min_b = 1
			b_incr = 3
			max_c = 150
			min_c = 50
			c_incr = 11
			max_d = 0.9
			min_d = 0.4
			d_incr = 6
			max_e = 0.9
			min_e = 0
			e_incr = 5.5
			start_solver_and_check_result(a, b, c, d, e, dataset, max_a, min_a, a_incr, max_b, min_b, b_incr, max_c, min_c, c_incr, max_d, min_d, d_incr, max_e, min_e, e_incr)

		if (dataset == 300): 
			a = 0.1;
			b = 1;
			c = 100;
			d = 0.7;
			e = 0;
			max_a = 0.3
			min_a = 0.1
			a_incr = 3
			max_b = 3
			min_b = 1
			b_incr = 3
			max_c = 150
			min_c = 50
			c_incr = 11
			max_d = 0.9
			min_d = 0.4
			d_incr = 6
			max_e = 0.9
			min_e = 0
			e_incr = 5.5
			start_solver_and_check_result(a, b, c, d, e, dataset, max_a, min_a, a_incr, max_b, min_b, b_incr, max_c, min_c, c_incr, max_d, min_d, d_incr, max_e, min_e, e_incr)

pool = mp.Pool(processes=int(num_cores / 2));
pool.map(launch_thread, datasets)
