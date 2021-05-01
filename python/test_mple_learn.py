import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
network = importr('network', on_conflict="warn",
                  robject_translations={
                      'as.tibble.network': 'as_tibble_network',
                      'as_tibble.network': 'as_tibble_network'
                  }
         )
base = importr('base')
ergm = importr('ergm')
readRDS = robjects.r['readRDS']
from rpy2.robjects import Formula
from mple_learn import stergmGraph




filename = '../data/edge_mutual__mutual.rds'
df = readRDS(filename)
df = np.array(df).astype(int)
n = df[0].shape[0]


form_terms = ['edges', 'mutual']
diss_terms = ['edges']


sim1 = stergmGraph(
    X = df,
    form_terms=form_terms,
    diss_terms=diss_terms)

yt0=df[0]
yt1=df[1]

lr_form, lr_diss = sim1.g_delta(yt0=yt0, yt1=yt1)

assert lr_form.shape == (n, n, len(form_terms))
assert lr_diss.shape == (n, n, len(diss_terms))

print("All test cases passed !!")