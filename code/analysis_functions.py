import numpy as np
import statsmodels.api as sm

# Functions for PyAtBSA analysis subprocess
def delta_snp_array(wtr, wta, mur, mua):
    #Calculate delta-SNP ratio
    return ((wtr)/(wtr+wta))-((mur)/(mur+mua))

def g_statistic_array(o1, o3, o2, o4):
    # Calculate G-statistic using numpy array input
    np.seterr(all='ignore')

    e1 = np.where(o1+o2+o3+o4!=0, (o1+o2)*(o1+o3)/(o1+o2+o3+o4), 0)
    e2 = np.where(o1+o2+o3+o4!=0, (o1+o2)*(o2+o4)/(o1+o2+o3+o4), 0)
    e3 = np.where(o1+o2+o3+o4!=0, (o3+o4)*(o1+o3)/(o1+o2+o3+o4), 0)
    e4 = np.where(o1+o2+o3+o4!=0, (o3+o4)*(o2+o4)/(o1+o2+o3+o4), 0)

    llr1 = np.where(o1/e1>0, 2*o1*np.log(o1/e1), 0.0)
    llr2 = np.where(o2/e2>0, 2*o2*np.log(o2/e2), 0.0)
    llr3 = np.where(o3/e3>0, 2*o3*np.log(o3/e3), 0.0)
    llr4 = np.where(o4/e4>0, 2*o4*np.log(o4/e4), 0.0)

    return np.where(e1*e2*e3*e4==0, 0.0, llr1+llr2+llr3+llr4)

def empircal_cutoff(posin, wt, mu):
    # Shuffle phenotype/genotype/position association 1000x
    # Establish variables
    lowess = sm.nonparametric.lowess
    lowess_span=0.3

    ## Establish lists
    smGstatAll = []
    smRatioAll = []
    RS_GAll = []
    smRS_G_yhatAll = []
    smRatio_yhatAll = []

    for i in range(1000):
        dfShPos = posin.sample(frac=1)
        dfShwt = wt.sample(frac=1)
        dfShmu = mu.sample(frac=1)

        smPos = dfShPos['pos'].to_numpy()
        sm_wt_ref = dfShwt['wt_ref'].to_numpy()
        sm_wt_alt = dfShwt['wt_alt'].to_numpy()
        sm_mu_ref = dfShmu['mu_ref'].to_numpy()
        sm_mu_alt = dfShmu['mu_alt'].to_numpy()
        
        #Calculate G-stats per iteration,collect results
        smGstat = g_statistic_array(sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt)
        smGstatAll.extend(smGstat)

        #Calculate snpRatio per iteration,collect results
        smRatio = delta_snp_array(sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt)
        smRatioAll.extend(smRatio)

        #Calculate ratio-scaled G-stat per iteration,collect results
        smRS_G = smRatio*smGstat
        RS_GAll.extend(smRS_G)

        #Lowess smooth per iteration, collect results
        smRS_G_yhatAll.extend(lowess(smRS_G,smPos, frac=lowess_span)[:,1])

    G_S_95p = np.percentile(smGstatAll, 95)
    RS_G_95p = np.percentile(RS_GAll, 95)
    RS_G_Y_99p = np.percentile(smRS_G_yhatAll, 99.99)

    return G_S_95p, RS_G_95p, RS_G_Y_99p