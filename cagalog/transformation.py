
from skbio.stats.composition import multiplicative_replacement
from skbio.stats.composition import clr
from sklearn import preprocessing
import numpy as np
import pandas as pd
import scipy.stats as ss

"""
Contains method for transforming data.

"""

def mult_replace(df):
    """
    wrapper for skbio's multiplicative multiplicative_replacement

    Parameters
    ----------
    df : DataFrame

    Returns
    -------
    df_mr : DataFrame
        modified via multiplicative replacement

    Notes
    -----
    Replaces zeros with the minimum non zero value in the entire
    matrix. Use multiplicaive replacement to ensure rows
    sum close to 1.

    """
    assert(isinstance(df, pd.DataFrame))
    nzra = np.min(df.values.flatten()[df.values.flatten() > 0])
    half_nzra = nzra/2
    # multiplicative replacement adds small value to non-zero entries while maintaining row sums equal to 1
    df_mr = pd.DataFrame(multiplicative_replacement(df, delta = half_nzra))
    assert(np.all(df_mr.values > 0))
    return df_mr

def box_cox(df):
    """
    Wrapper for sklearn's preprocessing.PowerTransformer (Box-Cox), which is
    only appropriate for strictly positive valued matrices. Zeros are replaced
    by this wrapper.

    Parameters
    ----------
    df : DataFrame


    Returns
    -------
    DataFrame
        Box-Cox transformed and
        modified via multiplicative replacement if zeros are present

    """
    assert(isinstance(df, pd.DataFrame))
    if not df.min().min() > 0:
        df = mult_replace(df)

    pt = preprocessing.PowerTransformer(method='box-cox',
                                        standardize=False)
    return pd.DataFrame(pt.fit_transform(df))

def yeo_johnson(df):
    """
    Wrapper for sklearn's preprocessing.PowerTransformer (Yeo-Johnson Option)
    which can handle negative values


    Parameters
    ----------
    df : DataFrame


    Returns
    -------
    DataFrame
        Yeo-Johnson transformed
    """
    assert(isinstance(df, pd.DataFrame))
    pt = preprocessing.PowerTransformer(method='yeo-johnson',
                                        standardize=False)
    return pd.DataFrame(pt.fit_transform(df))

def quantile_norm(df):
    """
    Wrapper for sklearn's preprocessing.QuantileTransformer.

    Parameters
    ----------
    df : DataFrame


    Returns
    -------
    DataFrame
        QuantileTransformer. transformed

    Notes
    -----
    Outer bounds are very low probability regions of the normal distribution so
    min and max are approximately -5 and +5 standard deviations away from mean,
    which limits the utility of this transform.
    """
    assert(isinstance(df, pd.DataFrame))
    qt = preprocessing.QuantileTransformer(output_distribution='normal',
                                           random_state=0)
    return pd.DataFrame(qt.fit_transform(df))

def rank_inv(df, stochastic = False):
    assert(isinstance(df, pd.DataFrame))
    """
    Wrapper for custom rank-based inverse normal transfromation that does put min and max
    so far away from the mean as preprocessing.QuantileTransformer.

    Parameters
    ----------
    df : DataFrames

    stochastic : boolean
        Set to false, otherwise ties and in particular zero values would be randomly ranked

    Notes
    -----
    This was written Ed Mountjoy (Statistical Geneticist at Open Targets,
    Wellcome Sanger Institute)
    code copied from:
    https://github.com/edm1/rank-based-INT (85cb37bb8e0d9e71bb9e8f801fd7369995b8aee7)

    """

    r = df.apply(lambda v : rank_INT(v, stochastic = stochastic), axis = 0)
    return pd.DataFrame(r)

def clr_with_mult_rep(df):
    assert(isinstance(df, pd.DataFrame))

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

    Notes
    -----
    This was written Ed Mountjoy (Statistical Geneticist at Open Targets,
    Wellcome Sanger Institute)
    code copied from:
    https://github.com/edm1/rank-based-INT (85cb37bb8e0d9e71bb9e8f801fd7369995b8aee7)

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
