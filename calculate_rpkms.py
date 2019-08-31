# Reads Per Kilobase Million as an alternative metric for using the CAGs data.
# apply to the commandline first argument the file path to hdf5 file

import sys
import pandas as pd
from cagalog.Cagalog import Cagalog
import cagalog.rpkm as rpkm #  see this module for rpkm functions.
import pickle

def main(fn):
    # Loading hdf5 file which contains CAG and metaphlan level data
    cg = Cagalog(hdf5_fp = fn)
    cg.read_cags()
    # gets matrix of PROP values (row : sample, columns: CAGs)
    cags = cg._pivot_cags()

    # Reads Per Kilobase Million DataFrames will be appended to d
    d  = [] # Contains rpkm by CAG
    #dg = [] # Contains rpkm by gene
    #da = [] # contains rpkm by allele

    # Loop through samples, and create a DataFame with rpk, rpkm. We also wanted
    # to compare this result with the PROP value calculated for CAGs.
    # Also for understanding the discrepancy between prop and rpkm,
    # we join the number of allele and the number of genes per CAG group to the
    # rpkm DataFrame.

    for sample in cags.index:
        print("Calculating RPKM for Sample {}".format(sample))
        r = rpkm.get_reads_per_sample(sample, fn)
        cag_rpkm = rpkm.calculate_rpkm_per_sample(r = r, sample_id = sample)
        gene_rpkm   = rpkm.calculate_gene_rpkm_per_sample(r = r, sample_id = sample)
        allele_rpkm = rpkm.calculate_allele_rpkm_per_sample(r = r, sample_id = sample)

        cag_rpkm_prop = pd.merge(cag_rpkm, cags.loc[sample,], how = 'left', left_on = 'group', right_index = True )
        x = gene_rpkm.copy().reset_index()
        y = allele_rpkm.copy().reset_index()
        genes_per_cag   = x.groupby(['group'])['gene'].count() # counts number of genes associated with each CAG group
        alleles_per_cag = y.groupby(['group'])['allele'].count() # counts number of alleles associated with each CAG group

        cag_rpkm_prop = pd.merge(cag_rpkm_prop,
                                 alleles_per_cag,
                                 how = "left",
                                 left_index = True,
                                 right_index = True)

        cag_rpkm_prop = pd.merge(cag_rpkm_prop,
                                 genes_per_cag,
                                 how = "left",
                                 left_index = True,
                                 right_index = True)

        # rename so all DataFrames have identical column names
        cag_rpkm_prop.rename(columns={ cag_rpkm_prop.columns[5]: "prop" }, inplace = True)
        d.append(cag_rpkm_prop)
        #dg.append(gene_rpkm)
        #da.append(allele_rpkm)
    return(d)




if __name__ == "__main__":
    fn = sys.argv[1] # apply to the commandline the file path to hdf5 file
    d = main(fn)
    # convert a list of your dataframes into one long dataframe
    combined_df = pd.concat(d, sort = False, axis = 0)
    pickle.dump( combined_df , open( "combined_df.p", "wb" ) )

    # This code is for plotting an discrepancy:

    # Here we observed an interesting discrepancy between rpkm and CAG - PROP.
    # When there is only 1 allele per gene the relationship is relatively linear;
    # however, when the allele per gene ratio increases PROP systematically
    # overestimates RPKM.
    d[0]['prop10k'] = 10000*d[0].iloc[:,5]
    d[0]['allele_per_gene'] = d[0]['allele'] /d[0]['gene']
    d[0].plot(kind = "scatter", x = 'rpk', y = 'prop10k', c = 'allele_per_gene', colormap = 'viridis' )
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages('rpk_vs_prop_for_CAGs.pdf') as pdf:
        d[0].plot(kind = "scatter", x = 'rpk', y = 'prop10k', c = 'allele_per_gene', colormap = 'viridis' )
        pdf.savefig()
