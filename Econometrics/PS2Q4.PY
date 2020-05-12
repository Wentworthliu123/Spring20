from numpy.random import multivariate_normal
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
import statsmodels.api as sm
from scipy.stats import percentileofscore
from scipy import stats
import seaborn as sns
import pandas as pd
from linearmodels.iv import IV2SLS


sim_N = 1000
    
N = 100
var_u = 1
var_v = 1
cov = 0.9
cov_matrix = [[var_u, cov],
             [cov, var_v]]
# DGP
beta = 1
pi_range = [2, 1, 0.25, 0]
beta_df = pd.DataFrame(columns=pi_range, index=['rejection ratio'])
for pi in pi_range:
# for pi in [2]:
    pi_result = []
    beta_result = []
    beta_t_result = []
    median = 0
    mean = 0
    for _i in range(sim_N):

        data = multivariate_normal([0,0], cov_matrix, size = N)
        u = data[:,0]
        v = data[:,1]
        z = np.random.randn(N)

        x = z*pi + v
        y = x*beta + u
        data=pd.DataFrame()
        data['x']=x
        data['y']=y
        data['z']=z
        model = sm.OLS(x, z).fit()
        pi_sim = model.params[0]
        mod = IV2SLS.from_formula('y ~ [ x ~ z]', data)
        beta_t = (mod.fit().params[0]-1)/mod.fit().std_errors[0]
        pi_result.append(model.tvalues[0])
        if np.abs(model.fvalue)>=10:
            beta_t_result.append(beta_t)
            beta_result.append(mod.fit().params[0])
        
    significant_pi = 100 - (percentileofscore(pi_result, 1.96) - percentileofscore(pi_result, -1.96))
    significant_beta = np.sum([1 for x in beta_t_result if np.abs(x)>=1.96])/len(beta_t_result)
    beta_df[pi] = significant_beta
    beta_df.loc['median',pi]=np.median(beta_result)
    beta_df.loc['mean',pi]=np.mean(beta_result)
    beta_df.loc['min',pi] = stats.describe(beta_result)[1][0]
    beta_df.loc['max',pi]=stats.describe(beta_result)[1][1]
    beta_df.loc['variance',pi]=stats.describe(beta_result)[3]
    beta_df.loc['skewness',pi]=stats.describe(beta_result)[4]
    beta_df.loc['kurtosis',pi]=stats.describe(beta_result)[5]+3
    

    print('For pi={}, {}% reject pi=0 at 5% level'.format(pi, round(significant_pi,2)))
    
#     pvalue = np.sum([1 for x in pvalue_result if x<0.05])
#     print('For pi={}, {}% reject significant_beta_t at 5% level'.format(pi, round(significant_beta_t,2)))
    print('{} reject beta=1 at 5% level, simulation median={}, mean={}\n'.format(round(significant_beta,2), round(median+1,4), round(mean+1,4)))

    ax = sns.distplot(beta_result,
                  bins=50,
                  kde=True,
                  hist_kws={"linewidth": 15,'alpha':1})
    plt.legend([2, 1, 0.25, 0])
    ax.set(xlabel='beta_hat distribution for each case', ylabel='Frequency')
    print(stats.describe(beta_result))
# plt.savefig("beta_distribution.png")  
beta_df