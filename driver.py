import sys
import os
import json
import datetime
import preprocessing
import pickle
import numpy as np
import time

t1 = datetime.datetime.now()
rdir = '/stor/tyche/hydro/private/lpt14/projects/routing_clustering/Fall2023/'
ncores = 16

# Read in the preprocessing metadata
metadata_file = '%s/metadata_datapreprocess.json' % rdir
metadata = json.load(open(metadata_file))

#Create the domain decomposition
preprocessing.domain_decomposition(metadata)

#Create the databases
os.system('mpirun -n %d python driver_core.py %s' % (ncores,metadata_file)) 
