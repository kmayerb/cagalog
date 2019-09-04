# cagalog - tools for working with the CAGalog

The CAGalog lives in hdf5.

The **Cagalog** class has some handy functions for accessing and transforming the data

* Centered log ratio transform is performed on full matrix with zeros replaced with
with half the lowest non-zero value across all samples.

* Prior to CLR transform - multiplicative replacement adds to the non-zero entries
while maintaining row sums equal to 1. (implemented with skbio.stats.composition)

* RPKM - calculate Reads Per Kb Million (i.e., the number of reads aligning to
a CAG, then adjusted for CAG size (kb) and divided the number total mapped
reads (in Millions) per sample)  

```python
from cagalog.Cagalog import Cagalog
import cagalog.rpkm as rpkm

fn = "path/filename.hdf5"
cg = Cagalog(hdf5_fp = fn)

cg.read_cags()
cg.read_metaphlan()

# get samples by cags matrix
cg._pivot_cags()

# get samples by taxa matrix at desired level of taxonomic resolution
cg._pivot_metaphlan("pyhlum")
cg._pivot_metaphlan("species")

# get CLR tranformed matrix
cg.clr_cags()
cg.clr_metaphlan("species")

# get  Reads per Kb Million (rpkm) by CAG, gene, or allele
r           = rpkm.get_reads_per_sample(sample, fn)
cag_rpkm    = rpkm.calculate_rpkm_per_sample(r = r, sample_id = sample)
gene_rpkm   = rpkm.calculate_gene_rpkm_per_sample(r = r, sample_id = sample)
allele_rpkm = rpkm.calculate_allele_rpkm_per_sample(r = r, sample_id = sample)
```

## Transformations

```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from cagalog.Cagalog import Cagalog
taxa = "family"
fn = 'path/filename.hdf5'
cg = Cagalog(fn)                                  # 1
cg.read_metaphlan()                               # 2
taxa_raw = cg._pivot_metaphlan(taxa)    # 3                   
ind_col = cg._filter_prevalence(taxa_raw, min_prev = .5).columns    # 4
taxa_clr = cg.clr_metaphlan(taxa)[ind_col]        # 5
taxa_raw = taxa_raw[ind_col]
bc = cg.transform(df = taxa_raw,  method = "box_cox") # 6
qn = cg.transform(df = taxa_raw , method = "quantile_norm")
ri = cg.transform(df = taxa_raw,  method = "rank_inv")
cl = cg.transform(df = taxa_raw , method = "clr")
dfs = [taxa_raw,cl,bc,qn,ri]
dfs[0].shape
```

```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from cagalog.Cagalog import Cagalog
from cagalog.Cagalog import RPKM
rpkm = RPKM(hdf5_fp = "cags_rpkm.hdf5")
str(rpkm)
rpkm.read_rpkm()
rpkm.pivot_rpkm()
fi = rpkm_filt = rpkm._filter_prevalence(min_prev= .9)
yj = rpkm.transform(method = "yeo_johnson")
bc = rpkm.transform(method = "box_cox")
qn = rpkm.transform(method = "quantile_norm")
ri = rpkm.transform(method = "rank_inv")
cl = rpkm.transform(method = "clr")
dfs = [fi.iloc[:, 1:10],
       cl.iloc[:, 1:10],
       yj.iloc[:, 1:10],
       qn.iloc[:, 1:10],
       ri.iloc[:, 1:10]]
cg.compare_dists(dfs)
```
