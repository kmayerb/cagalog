
from skbio.stats.composition import multiplicative_replacement
from skbio.stats.composition import clr
from sklearn import preprocessing
import numpy as np
import pandas as pd
import scipy.stats as ss

def mult_replace(df):
    """
    replace zeros with the minimum non zero value in the entire
    matrix. Use multiplicaive replacement to ensure rows
    sum close to 1.s
    """
    nzra = np.min(df.values.flatten()[df.values.flatten() > 0])
    half_nzra = nzra/2
    # multiplicative replacement adds small value to non-zero entries while maintaining row sums equal to 1
    df_mr = pd.DataFrame(multiplicative_replacement(df, delta = half_nzra))

    return(df_mr)

def box_cox(df):
    if not df.min().min() > 0:
        df = mult_replace(df)

    pt = preprocessing.PowerTransformer(method='box-cox',
                                        standardize=False)
    return pt.fit_transform(df)

def yeo_johnson(df):
    pt = preprocessing.PowerTransformer(method='yeo-johnson',
                                        standardize=False)
    return pt.fit_transform(df)

def quantile_norm(df):
    qt = preprocessing.QuantileTransformer(output_distribution='normal',
                                           random_state=0)
    return qt.fit_transform(df)

def rank_inv(df, stochastic = False):
    """
    Custom rank-based inverse normal transfromation that does put min and max
    so far away from the mean as preprocessing.QuantileTransformer
    """

    return df.apply(lambda v : rank_INT(v, stochastic = stochastic), axis = 0)

def clr_with_mult_rep(df):
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





def rank_INT(series, c=3.0/8, stochastic=True):
    """ Perform rank-based inverse normal transformation on pandas series.
        If stochastic is True ties are given rank randomly, otherwise ties will
        share the same value. NaN values are ignored.

        Args:
            param1 (pandas.Series):   Series of values to transform
            param2 (Optional[float]): Constant parameter (Bloms constant)
            param3 (Optional[bool]):  Whether to randomize rank of ties

        Returns:
            pandas.Series
    """

    # Check input
    assert(isinstance(series, pd.Series))
    assert(isinstance(c, float))
    assert(isinstance(stochastic, bool))

    # Set seed
    np.random.seed(123)

    # Take original series indexes
    orig_idx = series.index

    # Drop NaNs
    series = series.loc[~pd.isnull(series)]

    # Get ranks
    if stochastic == True:
        # Shuffle by index
        series = series.loc[np.random.permutation(series.index)]
        # Get rank, ties are determined by their position in the series (hence
        # why we randomised the series)
        rank = ss.rankdata(series, method="ordinal")
    else:
        # Get rank, ties are averaged
        rank = ss.rankdata(series, method="average")

    # Convert numpy array back to series
    rank = pd.Series(rank, index=series.index)

    # Convert rank to normal distribution
    transformed = rank.apply(rank_to_normal, c=c, n=len(rank))

    return transformed[orig_idx]

def rank_to_normal(rank, c, n):
    # Standard quantile function
    x = (rank - c) / (n - 2*c + 1)
    return ss.norm.ppf(x)
