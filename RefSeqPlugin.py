import PyPluMA
import PyIO

import sys
import os
import pandas as pd
import re
import numpy as np
import scipy as sp
from collections import defaultdict
from itertools import islice
os.environ['KMP_DUPLICATE_LIB_OK']='True'


class RefSeqPlugin:
 def input(self, inputfile):
  self.parameters = PyIO.readParameters(inputfile)
 def run(self):
  pass
 def output(self, outputfile):
  motifs = pd.read_csv(PyPluMA.prefix()+"/"+self.parameters["motifmap"], sep='\t', header=0)
  motifs.set_index('RBP', inplace=True, drop=False)
 
  rmap = defaultdict(dict)
  motif=""
  refs = {}
  with open(PyPluMA.prefix()+"/"+self.parameters["motifs"], "rt") as f:
    for l in f:
        if l.startswith('>'):
            l=re.sub('\s+$','',l)
            motif = l[1:]
            continue
        a = l.split('\t')
        rmap[motif][a[0]]=1
        refs[a[0]] = 1


  # In[13]:
  RDfile = PyPluMA.prefix()+"/"+self.parameters["RDfile"]#'input/high-vs-low_metastatic_lines_GSE59857_logFC_refseq.txt'

  exp = pd.read_csv(RDfile, sep='\t', header=0, index_col=0)


  # In[14]:


  mat = pd.DataFrame(0, index=list(set(exp.index) & set(refs.keys())), columns=motifs['RBP'])


  # In[15]:


  for rbp in mat.columns:
    m = motifs.loc[rbp,'motif']
    #sys.stdout.flush()
    ref = rmap[m].keys()
    ref = list(set(rmap[m].keys()) & set(mat.index))
    mat.loc[ref,rbp] = 1


  # In[16]:


  mat.to_csv(outputfile, sep='\t', index=True, index_label='RefSeq')

