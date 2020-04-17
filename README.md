# NEW: Volcanic.py

This repo is for developing tools for opening and interpretting the hdf5 output of [Golob-Minot geneshot](https://github.com/Golob-Minot/geneshot). The output of this NF workflow are documented in this [WIKI](https://github.com/Golob-Minot/geneshot/wiki/Getting-Started). 

## Curently 

cagalog.volcanic (volcanic.py) has tools for summarizing number of CAGs, Rank-Abundance, N50, and for making volcano plots from corncob output

Rvolcanic.R is R code for making more beautiful volcano plots

## Example Using IPython
```python
import os
from cagalog.volcanic import Seismic, StratoVolcano, Magma

fn_hdf5   = '/Volumes/LaCie/Users/kmayerbl/gscf/geneshot_cf_allfiles.results.hdf5'
fn_crncob = '/Volumes/LaCie/Users/kmayerbl/gscf/stats/corncob.results.csv'

assert os.path.isfile(fn_hdf5)
assert os.path.isfile(fn_crncob)

s = Seismic(fn_hdf5)
df = s._cag_size()

mg = Magma(fn_crncob)
mg.prep_volcano(var = 'mu.low_length', trim = 5000)
df_low = mg.magma['mu.low_length']['volcano']

cags = df_low.cag.to_list()
print("TOP HITS BY P-VALUE")
print(df_low.head(20))

sv = StratoVolcano(fn_hdf5)
# To get gene and tax_id for top hits
cag_details_df = sv._lookup_cag_list(cags[0:10])
print(cag_details_df )

# To investigate indivdual CAGs. E.g., the hits table anticorrelated with
negative_cags = df_low[df_low['est'] < 0]
sv._lookup_cag(9503) # Bifidobacterium longum
sv._lookup_cag(7992) # Enterobacteriaceae

# Or lookup the first 10 at once
cag_details_df = sv._lookup_cag_list(negative_cags.cag.head(10).to_list())
print(cag_details_df)

```

```
TRANQUILO, LOADING LOOKUP TABLES INTO MEMORY
NUMBER_OF_CAGS:            738631
NUMBER_OF_GENES:           4746411
NUMBER_OF_CAGS_W_ONE_GENE: 608999 (80.0)%
n50 (50% genes are on CAGs with >= 1001 genes)

RETURNING DataFrame:
CONSIDER RANK ABUNDANCE PLOT: df.plot.scatter(x = 'log10rank', y = 'log10size')
Please Be Tranquilo As We Measure The Input File Length
corncob.results.csv contains 25035810 lines
Progress[get variables]: [##########] 100% Done...Muito obrigado
VARIABLES DETECTED:
        mu.(Intercept)
        mu.time_pointm10
        mu.time_pointm12
        mu.time_pointm2
        mu.time_pointm3
        mu.time_pointm4
        mu.time_pointm5
        mu.time_pointm6
        mu.time_pointm8
        mu.time_pointm9
        mu.cf_statusControl
        mu.low_length
Progress[get coefs: mu.low_length]: [##########] 100% Done...Muito obrigado
Progress[get p_values: mu.low_length]: [##########] 100% Done...Muito obrigado
TOP HITS BY P-VALUE
           cag         pv       est
275772  316545  13.994784  0.793159
504493  579793  12.195587  1.383955
25482    29122  12.188693  0.714850
57578    65805  11.925344  0.693333
396343  455312  11.866951  1.320289
74249    84859  11.289885  0.760020
628249  722166  11.284005  1.379108
62213    71102  10.804078  0.747341
48239    55131  10.452281  0.646368
368251  422983  10.361562  1.017406
8088      9244  10.317737  0.595562
637550  732871  10.304813  0.751999
8315      9503  10.293219 -1.051609
7484      8554  10.244988  0.617622
6992      7992  10.232390 -0.864000
57706    65951  10.217337  0.656462
103215  118015  10.200530  0.736304
448645  515510  10.120851  0.836006
31565    36075  10.108978  0.724548
62946    71942  10.089501  0.669708
TRANQUILO, LOADING LOOKUP TABLES INTO MEMORY
('0-Name Not Found', 1)
('0-Name Not Found', 1)
('Enterococcus hirae', 8);('0-Name Not Found', 2)
('Lachnospiraceae', 4)
('Bifidobacterium breve', 1)
('Peptoniphilus harei', 2);('Peptoniphilus', 1)
('Enterococcus', 1)
('Fusobacterium', 2);('0-Name Not Found', 1)
('Klebsiella', 4);('0-Name Not Found', 1)
('Bifidobacterium longum', 1)
```


## Rank Abundance
![cag_rank_abundance](https://user-images.githubusercontent.com/46639063/79419208-3a073300-7f6b-11ea-98b0-cc84516206de.png)

## Show RVolcanic in Action


# DEPRECIATED: cagalog - tools for working with the CAGalog



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
