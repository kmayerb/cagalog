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
                if cnt % 100000 == 0:
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
                if cnt % 100000 == 0:
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
    def __init__(self, filename:str = None):
        if filename is None:
            self.filename = '/Volumes/LaCie/Users/kmayerbl/gscf/geneshot_cf_allfiles.result.hdf5'
        else:
            self.filename = filename

        print("TRANQUILO, LOADING BIG STUFF FOR YOU")
        self._lock_and_load()

    def _lock_and_load(self):
        """
        ONLY CALL THIS IF YOU NEED TO DO LOOKUPS. 
        """
        self.cag_to_gene    = pd.read_hdf(self.filename, "/annot/gene/cag")
        self.gene_to_taxid  = pd.read_hdf(self.filename,"/annot/gene/tax")
        self.taxid_to_annot = pd.read_hdf(self.filename, "/ref/taxonomy")
    
    def _lookup_cag(self,cag):
        df = self.cag_to_gene[self.cag_to_gene["CAG"] == cag].copy() 
        cag_genes = df.merge(self.gene_to_taxid, how = "left", left_on = "gene", right_on = "gene")
        return cag_genes
        #cags_genes_taxids = cag_genes.merge(self.taxid_to_annot, how = "left", left_on = "tax_id", right_on= "tax_id")
        #print(cags_genes_taxids.sample(10))
        #print(collections.Counter(cags_genes_taxids['tax_id']))
 



if __name__ == "__main__":


    import os
    #from cagalog import volcanic
    fn_hdf5 = '/Volumes/LaCie/Users/kmayerbl/gscf/geneshot_cf_allfiles.results.hdf5'
    assert os.path.isfile(fn_hdf5)
    # Get Summary of CAGs    
    # s = volcanic.Seismic(fn) 

    s = Seismic(fn)
    df = s._cag_size() 
    #df.plot.scatter(x = 'log10rank', y = 'log10size')

    fn = '/Volumes/LaCie/Users/kmayerbl/gscf/stats/corncob.results.csv'
    assert os.path.isfile(fn)

    #mg = volcanic.Magma(fn)
    mg = Magma(fn)
    mg.prep_volcano(var =  'mu.low_length', trim = 5000)
    #mg.prep_volcano(var =  'mu.cf_statusControl' , trim = 5000)
    df_low = mg.magma['mu.low_length']['volcano']
    #df_cf  = mg.magma['mu.cf_statusControl']['volcano']
    print(df_low.head(20))
    
    sv = StratoVolcano(fn_hdf5)
    cag = df_low.cag.iloc[0]
    sv._lookup_cag(cag)


    # IPYTHON MODE

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
    print(df_low.head(20))
    cags = df_low.cag.to_list()
    sv = StratoVolcano(fn_hdf5)
    for cag in cags[0:10]:
        cag = int(cag)
        print(cag)
        cag_df = sv._lookup_cag(cag)
        print(cag_df)


    #df.plot.scatter(x = 'log10rank', y = 'log10size')
    
    # mg._get_estimates(var =  'mu.low_length' )
    # mg._get_p_values(var =  'mu.low_length' )
    # df = mg._frame(var =  'mu.low_length')
    # df.plot.scatter(x = 'est',y = 'pv') 
    # mg._get_estimates(var =  'mu.cf_statusControl' )
    # mg._get_p_values(var =  'mu.cf_statusControl' )
    # df2 = mg._frame(var =  'mu.cf_statusControl')
    # df2.plot.scatter(x = 'est',y = 'pv') 

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


    

