import pandas as pd
import numpy as np
import requests
from collections import Counter
from time import sleep

class CAGinfo:
    """
    class specifically designed for getting cag info

    Attributes
    ----------
    fn : string
        filename of hdf5 used for lookup


    Methods
    -------
    lookup_all_taxids
        given a cag, looks up all associated taxids, and calls rest.ensemle.org to get
        taxonomic reference names

    get_taxids
        internal - for a given cag, gets the associated NCBI taxids for all the associated alleles

    tax_info
        internal - downloads json corresponding with taxid from  rest.ensemble.org

    tax_name
        internal - extracts taxonomic name form json downloaded from rest.ensemble.org


    Example
    -------
    >>> fn = "filename.hdf5"
    >>> caginfo = CAGinfo(fn)
    >>> cags = [9384, 11045, 46287]
    >>> {x:caginfo.lookup_all_taxids(x) for x in cags}
    {9384: [('Collinsella', 2),
    ('Faecalibacterium', 1),
    ('Clostridium sp. OM07-10AC', 1)],
     11045: [('Faecalibacterium', 3)],
     46287: [('Firmicutes', 1), ('Clostridium sp. AM34-11AC', 1)]}

    """
    def __init__(self,fn):
        """
        Init includes merger of CAGs and TAXID tables on 'allele' to enable CAG to taxid lookup

        Parameters
        ----------
        fn : string
            path and filename to .hdf5 file containing "groups/CAGs" and "groups/NCBI_TAXID"

        """
        self.fn = fn
        self.cags_df = pd.read_hdf(fn, "groups/CAGs")
        self.taxid_df = pd.read_hdf(fn, "groups/NCBI_TAXID")
        self.cags_to_taxid_df = self.cags_df.merge(self.taxid_df, on = ["allele"], how = "left")



    def lookup_all_taxids(self, cag):
        """
        Given a cag, looks up all associated taxids, and calls rest.ensemle.org to get
        taxonomic reference names

        Parameters
        ----------

        cag : int or string coercible to integer
            represents the cag group number

        Returns
        -------
        a Counter Dictionary with the taxid and number of alleles associated with it in that cag

        """
        taxids = self.get_taxids(cag)
        taxids = [x for x in taxids if int(x) is not 0]
        d = {x:self.tax_name(x) for x in set(taxids)}
        l = [d[x] for x in taxids]
        return Counter(l).most_common()

    def get_taxids(self, cag):
        """
        for a given cag, gets the associated NCBI taxids for all the associated alleles

        Parameters
        ----------

        cag : int or string coercible to integer
            represents the cag group number

        Returns
        -------
        an array of taxids associated with that group
        """
        if isinstance(cag, str):
            cag = int(cag)
        return self.cags_to_taxid_df.loc[self.cags_to_taxid_df["group"] == cag]["taxid"].values

    @classmethod
    def tax_name(self, taxid):
        """
        extracts taxonomic name form json downloaded from rest.ensemble.org

        Parameters
        ----------
        taxid : int

        Returns
        -------
        taxonomic_name : string

        """
        taxonomic_name = self.tax_info(taxid).get("name", str(taxid))
        return taxonomic_name

    @classmethod
    def tax_info(self, taxid, retries=10):

        """
        downloads json corresponding with taxid from  rest.ensemble.org"

        Parameters
        ----------
        taxid : int

        retries : int

        Returns
        -------
        r : json
            dictionary like object

        """
        taxid = str(taxid)
        r = requests.get("https://rest.ensembl.org/taxonomy/id/{}?content-type=application/json".format(taxid))
        n = 0
        while r.status_code != 200:
            sleep(1)
            n += 1
            if n > retries:
                print("Could not fetch {}".format(taxid))
                return {}
            r = requests.get("https://rest.ensembl.org/taxonomy/id/{}?content-type=application/json".format(taxid))

        return r.json()
