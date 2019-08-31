import os
import numpy as np
import pandas as pd
from skbio.stats.composition import multiplicative_replacement
from skbio.stats.composition import clr

class Cagalog:
    """
    Class for accessing metaphlan and cags data from the hdf5 format produced
    by the SM Lab.

    Attributes
    ----------
    hdf5_fp : string
        complete file path to hdf5 file

    min_prev_default : float
        the minimum prevalence across samples for a feature to be considered
        after filtering
        

    Methods
    -------
    read_metaphlan : method
        reads long form metaphlan data.frame from hdf5 tree

    read_cags
        reads long form cags data.frame from hdf5 tree

    clr_metaphlan()
        pivot + clr tranform with multiplicative replacement on metaphlan

    clr_cags()
        pivot + clr tranform with multiplicative replacement on cags

    _filter_prevalence
        filter column variable by minimum prevalence

    _pivot_metaphlan : method
        converts long to wide form using self.metaphlan limits
        to a specific taxonomic_level
    _pivot_cags(self, df=None):
        convert long to wide form using self.cags_df

    clr_transform_metaphlan_via_mult_rep_method : method
        uses multiplicative replacement to replace zeros with half of
        the lowest non-zero relative abundance value. Then performs clr
        transformation.

    clr_transform_cags_via_mult_rep_method : method
        uses multiplicative replacement to replace zeros with half of
        the lowest non-zero relative abundance value. Then performs clr
        transformation.

    Example
    -------
    import numpy as np
    import pandas as pd
    from Cagalog import Cagalog
    cg = Cagalog(hdf5_fp = "path/filename.hdf5")
    cg.read_cags()
    cg.read_metaphlan()
    cg.clr_transform_metaphlan_via_mult_rep_method("phylum")
    cg.clr_transform_metaphlan_via_mult_rep_method("family")
    cg.clr_transform_metaphlan_via_mult_rep_method("genus")
    cg.clr_transform_metaphlan_via_mult_rep_method("species")
    cg.clr_transform_cags_via_mult_rep_method()
    cg.cags_dict['cags_wide_mr_clr_df']

    """
    def __init__(self, hdf5_fp):
        self.min_prev_default = .5
        self.hdf5_fp = hdf5_fp  # complete file path to hdf5 file
        assert(os.path.isfile(self.hdf5_fp))

        self.x = 1               # for testing

        self.cags_dict = {}      # will hold raw and clr transformed CAGs result
        self.metaphlan_dict = {} # will hold raw and clr transformed metaphlan result
        self.cag_grouping_df = None


    # READER FUNCTIONS
    def read_metaphlan(self):
        """
        Reads long form metaphlan data.frame from hdf5 tree /abund/metaphlan/table"

        Assigns
        -------
        self.metaphlan_df : DataFrame

        """
        df = pd.read_hdf(self.hdf5_fp, "/abund/metaphlan/table")
        self.metaphlan_df = df

    def read_cags(self):
        """
        Reads long form cags data.frame from hdf5 tree /abund/CAGs"

        Assigns
        -------
        self.cags_df : DataFrame
        """
        df = pd.read_hdf(self.hdf5_fp,"/abund/CAGs")
        self.cags_df = df

    def read_CAGs_annotation(self):
        """
        Reads long form cags data.frame from hdf5 tree /groups/CAGs"
        """
        self.cag_grouping_df = pd.read_hdf(
            self.hdf5_fp,
            "/groups/CAGs"
        )
        return pd.read_hdf(self.hdf5_fp, "/groups/CAGs")

    def read_alleles(self):
        """
        "/annot/alleles"
        """
        return pd.read_hdf(self.hdf5_fp, "/annot/alleles")

    def read_taxid(self):
        """
        "/groups/NCBI_TAXID"

        """
        return pd.read_hdf(self.hdf5_fp, "/groups/NCBI_TAXID")

    # PIVOT FUNCTIONS
    def _pivot_metaphlan(self, taxonomic_level, df=None):
        """
        convert long to wide form using self.metaphlan limits
        to a specific taxonomic_level

        Returns
        -------
        mp_wide_taxa : DataFrame
        """
        if df is None:
            df = self.metaphlan_df
        mp = df[df['rank']== taxonomic_level]
        mp_wide_taxa = mp.pivot(index='sample',
                                columns='org_name',
                                values='values_block_0').fillna(0)
        return(mp_wide_taxa)

    def _pivot_cags(self, df=None):
        """
        convert long to wide form using self.cags_df

        Returns
        -------
        cags_wide : DataFrame
        """
        if df is None:
            df = self.cags_df
        cags_wide = df.pivot(index='sample',
                             columns='group',
                             values='prop').fillna(0)
        return(cags_wide)

    def _filter_prevalence(self, df, min_prev = None):
        """
        filter's column variable by min prevalence. i.e., X% percent of samples
        must have this variable wiht a value above 0 (or matrix min value)

        Parameters
        ----------
        df : DataFrame

        min_prev : float
            between 0 and 1

        Returns
        -------
        DataFrame
        """
        if min_prev is None:
            min_prev = self.min_prev_default
        min_val = df.min().min()
        to_keep = (df > min_val).mean() >= min_prev
        print("Features passing a {} prevalence threshold: {:,} / {:,}".format(
            min_prev ,
            to_keep.sum(),
            to_keep.shape[0]
        ))
        return df.loc[:, to_keep]


    def clr_metaphlan(self, taxonomic_level = "phylum"):
        """
        pivot + clr tranform with multiplicative replacement on metaphlan

        Parameters
        ----------
        taxonomic_level : string

        Returns
        -------
        DataFrame

        """
        return self._clr_transform_via_mult_rep_method(self._pivot_metaphlan("phylum"))

    def clr_cags(self):
        """
        pivot + clr tranform with multiplicative replacement on cags PROP variable

        Parameters
        ----------
        taxonomic_level : string

        Returns
        -------
        DataFrame

        """
        return self._clr_transform_via_mult_rep_method(self._pivot_cags())

    @classmethod
    def _clr_transform_via_mult_rep_method(self, df):
        nzra = np.min(df.values.flatten()[df.values.flatten() > 0])
        half_nzra = nzra/2
        # multiplicative replacement adds small value to non-zero entries while maintaining row sums equal to 1
        df_mr = multiplicative_replacement(df, delta = half_nzra)
        # clr transform
        mr_clr = clr(df_mr)
        # clr transform array to data.frame with index and column matching mp_wide_taxa
        mr_clr_df         = pd.DataFrame(mr_clr)
        mr_clr_df.columns = df.columns
        mr_clr_df.index   = df.index
        return mr_clr_df

    # clr transformation functions
    def clr_transform_metaphlan_via_mult_rep_method(self, taxonomic_level = "phylum"):
        """
        NOT GENERALIZABLE - DELETE

        uses multiplicative replacement to replace zeros with half of
        the lowest non-zero relative abundance value. Then performs clr
        transformation.

        Arguments
        ---------
        taxonomic_level : string
            "phlyum" through "species"

        Assigns
        -------
        self.metaphlan_dict : dictionary
            dictionary keyed on taxa level with the following attributes:
                1. mp_wide_taxa_df - taxa level relative abundances
                2. mp_wide_taxa_mr_clr_df - taxa level clr transformed
                   abundances (uses multiplicative replacement)
                3. half_nzra - on-zero relative abundance (NZRA) used for Mult Rep step

        """
        mp_wide_taxa = self._pivot_metaphlan(taxonomic_level = taxonomic_level)
        # one solution is to use the lowest non-zero relative abundance (NZRA), or more typically NZRA/2
        nzra = np.min(mp_wide_taxa.values.flatten()[mp_wide_taxa.values.flatten() > 0])
        half_nzra = nzra/2
        # multiplicative replacement adds small value to non-zero entries while maintaining row sums equal to 1
        mp_wide_taxa_mr = multiplicative_replacement(mp_wide_taxa, delta = half_nzra)
        # clr transform
        mp_wide_taxa_mr_clr = clr(mp_wide_taxa_mr)
        # clr transform array to data.frame with index and column matching mp_wide_taxa
        mp_wide_taxa_mr_clr_df         = pd.DataFrame(mp_wide_taxa_mr_clr)
        mp_wide_taxa_mr_clr_df.columns = mp_wide_taxa.columns
        mp_wide_taxa_mr_clr_df.index   = mp_wide_taxa.index

        self.metaphlan_dict[taxonomic_level] = {
                "mp_wide_taxa_df"        : mp_wide_taxa,
                "mp_wide_taxa_mr_clr_df" : mp_wide_taxa_mr_clr_df,
                "half_nzra"              : half_nzra }
        return(mp_wide_taxa_mr_clr_df)


    def clr_transform_cags_via_mult_rep_method(self):

        """
        NOT GENERALIZABLE - DELETE
        uses multiplicative replacement to replace zeros with half of
        the lowest non-zero relative abundance value. Then performs clr
        transformation.

        Arguments
        ---------
        taxonomic_level : string
            "phlyum" through "species"

        Assigns
        -------
        self.cags_dict : dictionary
        dictionary keyed on 'cags' with the following attributes:
            1. cags_wide_df -  relative abundances
            2. cags_wide_mr_clr_df - clr transformed
               abundances (uses multiplicative replacement)
            3. half_nzra - on-zero relative abundance (NZRA) used for Mult Rep step
        """

        cag_wide = self._pivot_cags()
        # one solution is to use the lowest non-zero relative abundance (NZRA), or more typically NZRA/2
        nzra = np.min(cag_wide.values.flatten()[cag_wide.values.flatten() > 0])
        half_nzra = nzra/2
        # multiplicative replacement adds small value to non-zero entries while maintaining row sums equal to 1
        cag_wide_mr = multiplicative_replacement(cag_wide, delta = half_nzra)
        # clr transform
        cag_wide_mr_clr = clr(cag_wide_mr)
        # clr transform array to data.frame with index and column matching mp_wide_taxa
        cag_wide_mr_clr_df         = pd.DataFrame(cag_wide_mr_clr)
        cag_wide_mr_clr_df.columns = cag_wide.columns
        cag_wide_mr_clr_df.index   = cag_wide.index

        self.cags_dict["cags"] = {
                "cags_wide_df"        :   cag_wide,
                "cags_wide_mr_clr_df" :   cag_wide_mr_clr_df,
                "half_nzra"          :   half_nzra }
        return cag_wide_mr_clr_df

        def fetch_metaphlan_result(self, clr= True, taxonomic_level = "phylum"):
            """
            getter
            """
            if clr:
                key = 'mp_wide_taxa_mr_clr_df'
            else:
                key = 'mp_wide_taxa_df'
            try:
                return(self.metaphlan_dict[taxonomic_level][key])
            except KeyError:
                print("NO METAPHLAN MATRIX CREATED SEE clr_transform_metaphlan_via_mult_rep_method()")




#from collections import defaultdict
#from functools import lru_cache
#import json
#import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt
#import seaborn as sns
#import scipy
#import sklearn
#from time import sleep
#from random import shuffle
#from sklearn.decomposition import PCA
#import statsmodels.api as sm
#import os
#from statsmodels.stats.multitest import multipletests
#from scipy.spatial.distance import pdist, squareform
