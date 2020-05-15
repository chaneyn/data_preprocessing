import sys
import os
import numpy as np
import json
import pickle
import preprocessing
import gc
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print('here',flush=True)
#name = MPI.Get_processor_name()
#rank = 0#12
#size = 1
#size = 29#16
#rank = 39#106
#size = 107
eares = 26.0#670.0**0.5 #ONLY APPLICABLE FOR THE GIVEN REGION. NEED TO GENERALIZE

#Get general info
pre_metadata_file = sys.argv[1]

#Read in the metadata
pre_metadata = json.load(open(pre_metadata_file))

#Read in the catchment summary database
pck_file = '%s/domain/domain_database.pck' % pre_metadata['output_data']
cdb = pickle.load(open(pck_file,'rb'))
crange = range(len(cdb))

for ic in crange[rank::size]:

 print("Rank:%d, Catchment:%s - Initializing" % (rank,ic),flush=True)

 cid = cdb[ic]['cid']
 #cid = 1

 #Define the catchment directory
 cdir = '%s/domain/%d' % (pre_metadata['output_data'],cid)

 # Update Metadata
 '''metadata['Preprocessing']['cid'] = str(cid)
 metadata['Preprocessing']['workspace']   = 'input_data/catchments/%d/workspace' % cid
 metadata['Preprocessing']['input_file']  = 'output_data/%d/input_file.nc' % cid
 metadata['HydroBlocks']['workspace']   = 'input_data/catchments/%d/workspace' % cid
 metadata['HydroBlocks']['input_file']  = 'output_data/%d/input_file.nc' % cid
 metadata['HydroBlocks']['output']['dir'] = 'output_data/%d' % cid
 metadata['HydroBlocks']['restart']['dir'] = 'output_data/%d/restart' % cid
 json.dump(metadata,open('metadata_%d.json' % cid ,'w'))'''
   
 #Prepare the workspace for the sub-domain
 os.system('rm -rf %s' % cdir)
 os.system('mkdir -p %s' % cdir)
 preprocessing.prepare_input_data(cdir,cdb[ic],pre_metadata,rank,ic,eares)
 #exit()
