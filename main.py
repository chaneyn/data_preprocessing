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

#Get general info
pre_metadata_file = sys.argv[1]

#Read in the metadata
pre_metadata = json.load(open(pre_metadata_file))

#Create the domain decomposition
if rank == 0:
  preprocessing.domain_decomposition(pre_metadata)
comm.Barrier()

#Correct and finalize the domain decomposition
preprocessing.correct_domain_decomposition(comm,pre_metadata)
comm.Barrier()

#Read in the catchment summary database
pck_file = '%s/cids/domain_database.pck' % pre_metadata['output_data']
cdb = pickle.load(open(pck_file,'rb'))
crange = range(len(cdb))

for ic in crange[rank::size]:

 print("Rank:%d, Catchment:%s - Initializing" % (rank,ic),flush=True)

 cid = cdb[ic]['cid']
 #cid = 1

 #Define the catchment directory
 cdir = '%s/cids/%d' % (pre_metadata['output_data'],cid)

 #Prepare the workspace for the sub-domain
 os.system('rm -rf %s' % cdir)
 os.system('mkdir -p %s' % cdir)
 preprocessing.prepare_input_data(cdir,cdb[ic],pre_metadata,rank,ic)
