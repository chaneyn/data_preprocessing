import sys
import os
import json
import datetime
import preprocessing
import pickle
import numpy as np
import time

#time.sleep(5000)

t1 = datetime.datetime.now()
#rdir = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/regions/LW_OK'
rdir = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/regions/LR_GA'
#rdir = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/regions/SGP'
#rdir = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/regions/yellowstone'
ncores = 9#256

# Read in the preprocessing metadata
metadata_file = '%s/metadata_data_preprocess.json' % rdir
metadata = json.load(open(metadata_file))

#Create the domain decomposition
preprocessing.domain_decomposition(metadata)

#Create the databases
#print('here1',flush=True)
#print('mpirun -n %d python driver_core.py %s' % (ncores,metadata_file))
os.system('mpirun -n %d python driver_core.py %s' % (ncores,metadata_file)) 
#os.system('python driver_core.py %s' % (metadata_file,)) 
