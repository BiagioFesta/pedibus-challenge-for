#!/bin/python
# Copyright 2017 <Biagio Festa>
import multiprocessing as mp
import os
import subprocess as sp

root_dir = "."
solver = "../for-ch_solver"
num_cores = mp.cpu_count()
list_dataset = [f for f in os.listdir(root_dir) if f.endswith(".dat")]

print("Num cores detected: " + str(num_cores))
print("Dataset files detected:", list_dataset)
print("Solver detected: " + solver)

def exec_dataset(file):
    param = "-t 3600"
    command = solver + " " + str(file) + " " + param
    print("Execution:", command)
    sp.call(command, shell = True)
    

pool = mp.Pool(processes=int(num_cores / 2));
pool.map(exec_dataset, list_dataset)
