"""
make volcano plots without reading entire file into memory
"""

import sys
import os 
import collections
import pandas as pd
import numpy as np


class Seismic():
    """
    Before making volcano plots, consider number of CAGs and distribution of sizes 

    Example
    -------
    In [1]:     from cagalog import volcanic 
       ...:     s = volcanic.Seismic() 
       ...:     df = s._cag_size() 
       ...:     df.plot.scatter(x = 'log10rank', y = 'log10size')    

    NUMBER_OF_CAGS:            738631
    NUMBER_OF_GENES:           4746411
    NUMBER_OF_CAGS_W_ONE_GENE: 608999 (80.0)%
    n50 (50% genes are on CAGs with >= 1001 genes)
    
    RETURNING DataFrame:
    CONSIDER RANK ABUNDANCE PLOT: df.plot.scatter(x = 'log10rank', y = 'log10size')
    Out[1]: <matplotlib.axes._subplots.AxesSubplot at 0x11c7ce588>
    """
    def __init__(self, filename:str = None):
        if filename is None:
            self.filename = '/Volumes/LaCie/Users/kmayerbl/gscf/geneshot_cf_allfiles.results.hdf5'
        else:
            self.filename = filename
    
    def _get_n50(self,x,n:int,f:float = .5):
        """
        x : Series 
            (must be sorted in descending order)
        n : int 
            number of total genes 
        f : float 
            .5 for N50
        """
        n #len(cag_to_gene)
        nx = n * f
        xcs = x.cumsum()
        xcsbool = xcs > nx
        i = xcsbool.idxmax()
        return x[i]

    def _cag_size(self, dest:str = "/annot/gene/cag"):
        """
        dest : string 
            destination in the hdf5 
        """
        cag_to_gene  = pd.read_hdf(self.filename, dest)
        x = cag_to_gene['CAG'].value_counts() 
        
        n_genes = np.sum(x)
        n_cags = len(x)
        n_singletons = np.sum(x == 1)
        p_singletons = 100*round(n_singletons/n_cags , 1)    
        x = x.sort_index()
        n50 = self._get_n50(x=x, n=n_genes, f = .5)
        sys.stdout.write(f'NUMBER_OF_CAGS:            {n_cags}\n')
        sys.stdout.write(f'NUMBER_OF_GENES:           {n_genes}\n')
        sys.stdout.write(f'NUMBER_OF_CAGS_W_ONE_GENE: {n_singletons} ({p_singletons})%\n')  
        sys.stdout.write(f'n50 (50% genes are on CAGs with >= {n50} genes)\n\n')
        df = pd.DataFrame({"cag":x.index, "size":x})
        df["log10size"] = df['size'].apply(np.log10)
        df["log10rank"] = df['cag'].add(1).apply(np.log10)
        sys.stdout.write('RETURNING DataFrame:\n')    
        sys.stdout.write("CONSIDER RANK ABUNDANCE PLOT: df.plot.scatter(x = 'log10rank', y = 'log10size')\n")
        return df



class Magma:
    """
    Magma makes volcanos
    """
    def __init__(self, filename:str, formula = None):
        self.filename  = filename
        self.formula   = formula 
        self.variables = None
        sys.stderr.write("Please Be Tranquilo As We Measure The Input File Length\n")
        self._get_nlines()
        sys.stderr.write(f"{os.path.basename(self.filename)} contains {self.nlines} lines\n")
        self.magma     = dict()
        self._get_variables()

    def prep_volcano(self,var, trim = 1000):
        self._get_estimates(var = var )
        self._get_p_values(var = var )
        self._frame_df(var = var , trim = trim) 

    def _get_nlines(self):
        cmd = f'wc -l {self.filename} > nlines.tmp'
        os.system(cmd)
        with open('nlines.tmp', 'r') as fh:
            n = fh.readlines()[0].strip().split(" ")[0]
            n = int(n)
        self.nlines = n

    def _get_variables(self):
        c = collections.Counter()
        self._update_progress(0, task_str = "get variables")
        cnt = 0
        with open(self.filename, "r") as fh:
            fh.readline()
            for line in fh:
                cnt += 1
                #if cnt % 100 == 0:
                    #self._update_progress(cnt/float(self.nlines), task_str = "get variables")
                parameter, param_type, value, cag = line.strip().split(",")
                if param_type == "estimate":
                    c.update({parameter:1})
                if int(cag) > 10: 
                    break
        self._update_progress(1,task_str = "get variables")
        self.variables = list(c.keys())  
        sys.stdout.write("VARIABLES DETECTED:\n\t")
        sys.stdout.write("\n\t".join(map(str, list(c.keys()))))
        return(self.variables)
    
    def _get_estimates(self, var=None):
        self._update_progress(0,task_str = f"get coefs: {var}")
        cnt = 0
        with open(self.filename, "r") as fh:
            d = collections.OrderedDict()
            for line in fh:
                cnt += 1
                if cnt % 1000000 == 0:
                    self._update_progress(cnt/float(self.nlines),task_str = f"get coefs: {var}")
                parameter, param_type, value, cag = line.strip().split(",")
                if param_type == "estimate" and parameter == var:
                    try:
                        d.update({cag:float(value)})
                    except ValueError:
                        d.update({cag: None})
        self._update_progress(1,task_str = f"get coefs: {var}")
        if var not in self.magma.keys():
            self.magma[var] = dict()
        self.magma[var]['estimates'] = d
           
    def _get_p_values(self, var):
        cnt = 0
        self._update_progress(0,task_str = f"get p_values: {var}")
        with open(self.filename, "r") as fh:
            d = collections.OrderedDict()
            for line in fh:
                cnt += 1
                if cnt % 1000000 == 0:
                    self._update_progress(cnt/float(self.nlines),task_str = f"get p_values: {var}")
                parameter, param_type, value, cag = line.strip().split(",")
                if param_type == "p_value" and parameter == var:
                    try:
                        d.update({cag:float(value)})
                    except ValueError:
                        d.update({cag: None})
        self._update_progress(1,task_str = f"get p_values: {var}")
        if var not in self.magma.keys():
            self.magma[var] = dict()
        self.magma[var]['p_values'] = d
   
    def _frame_df(self, var:str, sort:bool = True, trim:int = None):
        pv = self.magma[var]['p_values']  
        est = self.magma[var]['estimates']
        dfpv  = pd.DataFrame([(k,v) for k,v in pv.items()], columns = ['cag', 'pv']) 
        dfest = pd.DataFrame([(k,v) for k,v in est.items()], columns = ['cag', 'est'])
        df = pd.merge(dfpv, dfest, how = 'left', left_on = 'cag', right_on = 'cag') 
        if sort:
            df = df.sort_values(['pv']) 
            df['pv'] = df['pv'].apply(np.log10).multiply(-1).copy()   
        if trim is None:
            self.magma[var]['volcano'] = df.copy()
        else:
            self.magma[var]['volcano'] = df.iloc[0:trim,:]



    def _scan(self):
        with open(self.filename, "r") as fh:
            c = 0
            for line in fh:
                parameter, param_type, value, cag = line.strip().split(",") 
                if param_type == "estimate":
                    print(parameter, param_type, value, cag )
                
                    c += 1
                    if c > 5:
                        break 

    def _enumerate(self, type:str = "estimate"):
        pass
        
    def _identify_model(self):
        pass

    def _update_progress(self, progress, task_str =""):
        """
        https://stackoverflow.com/questions/3160699/python-progress-bar
        """
        barLength = 10 # Modify this to change the length of the progress bar
        status = ""
        if isinstance(progress, int):
            progress = float(progress)
        if not isinstance(progress, float):
            progress = 0
            status = "error: progress var must be float\r\n"
        if progress < 0:
            progress = 0
            status = "Halt...\r\n"
        if progress >= 1:
            progress = 1
            status = "Done...Muito obrigado\r\n"
        block = int(round(barLength*progress))
        progress = round(progress*100,2)
        text = "\rProgress[{0}]: [{1}] {2}% {3}".format(  task_str, "#"*block + "-"*(barLength-block), progress, status)
        sys.stdout.write(text)
        sys.stdout.flush()


class StratoVolcano():
    """
    Class for evaluating individual CAG layers discovered 
    from a volcano analysis.

    Example
    -------
    import os
    from cagalog.volcanic import Seismic, StratoVolcano, Magma 
    fn_hdf5   = '/Volumes/LaCie/Users/kmayerbl/gscf/geneshot_cf_allfiles.results.hdf5' 
    sv = StratoVolcano(fn_hdf5) 
    sv._lookup_cag(1) 
    """
    def __init__(self, filename:str = None):
        """
        filename : str
            full path to the GeneShot Results File
        """
        if filename is None:
            self.filename = '/Volumes/LaCie/Users/kmayerbl/gscf/geneshot_cf_allfiles.result.hdf5'
        else:
            self.filename = filename

        # <fn_taxid_dict> is the full path to lookup table mapping
        # taxid as sciname
        self.fn_taxid_dict = "/Users/kmayerbl/active/cagalog/taxdump/taxid_to_sciname.csv"

        # Initialize Empty Objects
        # <cag_to_gene> DataFrame (source: "/annot/gene/cag" )
        self.cag_to_gene = None
         # <gene_to_taxid_dict> dictionary (source: "/annot/gene/tax")
        self.gene_to_taxid_dict = None 
        # <taxid_to_sciname_dict> dictionary (source: loaed from csv e.g., self.fn_taxid_dict)
        self.taxid_to_sciname_dict = None

        # Load Objects
        sys.stderr.write("TRANQUILO, LOADING LOOKUP TABLES INTO MEMORY\n")
        self._load_cag_to_gene_df()
        self._load_gene_to_taxid_dict()      
        self._load_taxid_to_sciname_dict()
        
    def _load_gene_to_taxid_dict(self):
        gene_to_taxid  = pd.read_hdf(self.filename,"/annot/gene/tax")
        self.gene_to_taxid_dict = {k:v for k,v in zip(gene_to_taxid.gene, gene_to_taxid.tax_id)}
        self.gene_to_taxid = gene_to_taxid

    def _load_cag_to_gene_df(self):
        """
        LOAD TO MEMORY FOR FAST LOOKUPS
        """
        self.cag_to_gene    = pd.read_hdf(self.filename, "/annot/gene/cag")    

    def _load_taxid_to_sciname_dict(self, fn = None):
        if fn is None:
            fn = self.fn_taxid_dict 
        taxid_to_sciname_dict = dict()
        with open(fn, "r") as fh:
            for line in fh:
                try:
                    k,v = line.strip().split(",")
                except ValueError:
                    continue
                taxid_to_sciname_dict[k] = v
        self.taxid_to_sciname_dict = taxid_to_sciname_dict 

    def _safely_taxid_to_sciname_dict(self,k):
        try:
            return self.taxid_to_sciname_dict[str(k)]
        except KeyError:
            return f"{k}-Name Not Found"

    def _safely_gene_to_taxid_dict(self, k):
        try:
            return self.gene_to_taxid_dict[str(k)]
        except KeyError:
            return f"{k}-Gene Not Found"

    def _lookup_cag(self, cag, verbose = False):
        
        df = self.cag_to_gene[self.cag_to_gene["CAG"] == cag].copy() 

        df['tax_id'] = df['gene'].apply(lambda k: self._safely_gene_to_taxid_dict(k))
        
        # provide scientific names from a pre-loaded dictionary       
        df['sciname'] = df['tax_id'].\
                                apply(lambda k: self._safely_taxid_to_sciname_dict(k))
        if verbose == True:
            sys.stdout.write(";".join(map(str, collections.Counter(df['sciname']).most_common())) + "\n" )
       
        return df
         
    def _lookup_cag_taxonomy_string(self,cag:int)->str:
        df = self._lookup_cag(cag = cag, verbose =False)
        ts = ";".join(map(str, collections.Counter(df['sciname']).most_common() ))
        return ts 
    
    def _lookup_cag_genes_string(self,cag:int)->str:
        df = self._lookup_cag(cag = cag, verbose =False)
        ts = ";".join(map(str, df['gene']))
        return ts 

    def _lookup_cag_taxonomy_dict(self,cag:int)->dict:
        df = self._lookup_cag(cag = cag, verbose =False)
        return {cag : collections.Counter(df['sciname']).most_common()}

    def _lookup_cag_list(self, cags:list):
        """
        cags : list of ints
            list of cag numbers
        """
        cag_dfs = list()
        for cag in cags:
            if not isinstance(cag, int):
                cag = int(cag)
            cag_df = self._lookup_cag(cag)
            ts = self._lookup_cag_taxonomy_string(cag)
            gs = self._lookup_cag_genes_string(cag)
            most_common_ts = collections.Counter(cag_df.sciname).most_common()[0][0] 
            cag_df['ts_mode'] = most_common_ts 
            cag_df['ts'] = str(ts)
            cag_df['gs'] = str(gs)
            cag_df = cag_df.rename(columns = {"CAG":"cag"})
            cag_dfs.append(cag_df)
        final_cag_df = pd.concat(cag_dfs).reset_index(drop = True)
        return final_cag_df

    def _lookup_cag_list_summary(self, cags : list):
        df = self._lookup_cag_list(cags)
        df = df.groupby(["cag"]).first().reset_index()   
        df = df[['cag','ts_mode','ts','gs']].copy()
        return(df)


if __name__ == "__main__":

    import os
    from cagalog.volcanic import Seismic, StratoVolcano, Magma
    import pandas as pd

    fn_hdf5   = '/Volumes/LaCie/Users/kmayerbl/gscf/geneshot_cf_allfiles.results.hdf5'
    fn_crncob = '/Volumes/LaCie/Users/kmayerbl/gscf/stats/corncob.results.csv'

    assert os.path.isfile(fn_hdf5)
    assert os.path.isfile(fn_crncob)

    s = Seismic(fn_hdf5)
    df = s._cag_size()

    mg = Magma(fn_crncob)
    mg.prep_volcano(var = 'mu.low_length', trim = 5000)
    df_low= mg.magma['mu.low_length']['volcano']

    cags = df_low.cag.to_list()
    print("TOP HITS BY P-VALUE")
    print(df_low.head(20))

    sv = StratoVolcano(fn_hdf5)
    # To get gene and tax_id for top hits
    cag_details_df = sv._lookup_cag_list(cags[0:10])
    print(cag_details_df )

    # To investigate indivdual CAGs. E.g., the hits table anticorrelated with
    negative_cags = df_low[df_low['est'] < 0]
    positive_cags = df_low[df_low['est'] > 0]

    negative_cags = negative_cags.apply(pd.to_numeric)
    positive_cags  = positive_cags.apply(pd.to_numeric) 

    sv._lookup_cag(9503) # Bifidobacterium longum
    sv._lookup_cag(7992) # Enterobacteriaceae

    # Or lookup the first 10 at once
    cag_neg_details_df = sv._lookup_cag_list_summary(negative_cags.cag.head(20).to_list())
    cag_pos_details_df = sv._lookup_cag_list_summary(positive_cags.cag.head(20).to_list())
    negative_cags['cag'] = negative_cags.cag.astype(int).to_list()

    negative_cags.head(20).merge(cag_neg_details_df , how = "left", right_on = "cag", left_on= "cag").to_csv('low.weight.negative20.csv')
    positive_cags.head(20).merge(cag_pos_details_df , how = "left", right_on = "cag", left_on= "cag").to_csv('low.weight.positive20.csv')

    # # another easy way with pandas
    # fn = '/Volumes/LaCie/Users/kmayerbl/gscf/stats/corncob.results.csv'
    # import pandas as pd
    # import numpy as np
    # df = pd.read_csv(fn)
    # ind = df['parameter'] == 'mu.cf_statusControl'
    # dft = df[ind]
    # dftu = dft.set_index(['parameter', 'CAG',"type"]).unstack().reset_index()
    # dftu.columns = dftu.columns.map(''.join)
    # dftu['log10_pv'] = dftu['valuep_value'].apply(np.log10).multiply(-1)
    # dftu.plot.scatter(x = 'valueestimate', y = 'log10_pv', alpha = .1 )


    

