import sys
import numpy as np
import pandas as pd
from cagalog.Cagalog import Cagalog

def get_reads_per_sample(sample_id, local_hdf5):
    """
    Gets depth, nreads, and length of all alleles from a given sample

    Parameters
    ----------
    sample_id : string
    local_hdf5 : string

    Returns
    -------
    cag_df : pd.DataFrame

    # read in the abundance for one particular sample in this dataset
    # read in the grouping of alleles for all cags
    # add the abundance information to each of the alleles from this CAG
    # remove the alleles which weren't observed in this sample
    # return cag_df
    """
    sample_allele_abund = pd.read_hdf(local_hdf5,
                                      "/abund/alleles",
                                      where="sample == '{}'".format(sample_id)).set_index("id")
    cag_df = pd.read_hdf(local_hdf5, "/groups/CAGs")
    for k in ["depth","nreads", "length"]:
        cag_df[k] = cag_df["allele"].apply(sample_allele_abund[k].get)
    cag_df = cag_df.dropna()
    return cag_df

def calculate_rpkm_per_sample(r, sample_id):
    """
    Use allele information from 1 Sample to calculate RPKM (Reads per Kb Million)
    per CAG!

    Parameters
    ----------
    r : DataFrame
        containing allele information from 1 sample at a minimum ["group", "nreads", "length" ]
    sample_id : string

    Returns
    -------
    rs : DataFrame
        contains 'group' (CAG group), 'rpk' (reads per Kb), and 'rpkm' (reads per Kb Millon)

    Notes
    -----
    # 1. calculate total reads per sample
    # 2. divide total reads per sample by 1 Million
    # 3. sum reads by cag group, sum lengths by cag
    # 4. rpk =   #reads per cag  /  length of all alleles comprising cag
    # 5. rpkm =  rpk / totall_mapped_reads in millions

    """
    total_mapped_million = r['nreads'].sum() / 1000000
    rs = r.groupby(['group'])['nreads','length'].sum()
    rs['rpk'] = rs['nreads'] / rs['length']
    rs['rpkm'] = rs['rpk'] / total_mapped_million
    rs['sample'] = sample_id
    return(rs)

def calculate_gene_rpkm_per_sample(r, sample_id):
    """
    Use allele information from 1 Sample to calculate RPKM (Reads per Kb Million)
    per gene!

    Parameters
    ----------
    r : DataFrame
        containing allele information from 1 sample at a minimum ["group", "nreads", "length" ]
    sample_id : string

    Returns
    -------
    rs : DataFrame
        contains 'group' (CAG group), 'gene' (gene), 'rpk' (reads per Kb), and 'rpkm' (reads per Kb Millon)

    Notes
    -----
    # 1. calculate total reads per sample
    # 2. divide total reads per sample by 1 Million
    # 3. sum reads by cag group and gene sum lengths effectively by gene
    # 4. rpk =   #reads per gene /  length of all alleles comprising a gene
    # 5. rpkm =  rpk / totall_mapped_reads in millions

    """
    total_mapped_million = r['nreads'].sum() / 1000000
    rs = r.groupby(['gene',"group"])['nreads','length'].sum()
    rs['rpk'] = rs['nreads'] / rs['length']
    rs['rpkm'] = rs['rpk'] / total_mapped_million
    rs['sample'] = sample_id
    return(rs)

def calculate_allele_rpkm_per_sample(r, sample_id):
    """
    Use allele information from 1 Sample to calculate RPKM (Reads per Kb Million)
    per allele!
    """
    total_mapped_million = r['nreads'].sum() / 1000000
    rs = r.groupby(['allele','gene', "group"])['nreads','length'].sum()
    rs['rpk'] = rs['nreads'] / rs['length']
    rs['rpkm'] = rs['rpk'] / total_mapped_million
    rs['sample'] = sample_id
    return(rs)
