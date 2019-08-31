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
