import pandas as pd
import numpy as np
import scipy
import control
from statsmodels.tsa.vector_ar.var_model import VAR
import statsmodels.api as sm
from typing import Optional

class EsmssmStats:

    def __init__(self):
        pass
        
    
    def pval_to_asterix(self, p: float, multiple_comparison: Optional[float]=None) -> str:
        if multiple_comparison == None: multiple_comparison = 1
        if p <= 0.001/multiple_comparison: asterix = '***'
        elif p <= 0.01/multiple_comparison: asterix = '**'
        elif p <= 0.05/multiple_comparison: asterix = '*'
        else: asterix = '$^{ns}$'
        return asterix

    
    def replace_small_pvalues(self, df: pd.DataFrame) -> pd.DataFrame:
        sig = df['pvalue'] < 0.001
        df.loc[~sig, 'pvalue'] = df['pvalue'][~sig].round(3)
        df.loc[sig, 'pvalue'] = '$<$0.001'
        return df

    
    def test_normality(self, data: np.array) -> list:
        dim = np.shape(data)
        normality_assumption = []
        if len(dim) > 2:
            for i in range(dim[0]):
                for j in range(dim[1]):
                    normality_test = scipy.stats.shapiro(data[i,j,:])
                    normality_assumption.append(normality_test[1] > 0.05)
        elif len(dim) == 2:
            for i in range(dim[0]):
                normality_test = scipy.stats.shapiro(data[i,:])
                normality_assumption.append(normality_test[1] > 0.05)
        elif len(dim) == 1:
            normality_test = scipy.stats.shapiro(data)
            normality_assumption.append(normality_test[1] > 0.05)

        return normality_assumption
        

    def fit_simulate_var(self, data: np.array) -> tuple:
        T,D = np.shape(data)
        var = VAR(data)
        results = var.fit(1)
        forcast = results.forecast(data, steps=T) + np.random.multivariate_normal(np.zeros(D), np.cov(results.resid.T), T)
        return results, forcast
        

    def extract_parameter_from_mat(self, data: np.array) -> dict:
        params = dict()
        params['B'] = data['B'][0][0]
        params['A'] = data['A'][0][0] + data['W'][0][0]
        params['h'] = np.squeeze(data['h'][0][0])
        params['C'] = data['C'][0][0]
        params['Q'] = data['S'][0][0]
        params['R'] = data['G'][0][0]
        params['mu_0'] = np.squeeze(data['mu0'][0][0])
        params['sigma_0'] = data['S'][0][0]                         
        return params
        

    def simulate_KF(self, data: np.array, parameter: dict) -> np.array:
        options = dict()
        options['inputs'] = False
        params = self.extract_parameter_from_mat(parameter)
        z,x = self.sample_lds(len(data), params, options)
        return x
        

    def sample_lds(self, n_timesteps: int, params: dict, options: dict) -> tuple:
        """
        Generate samples from a Linear Dynamical System.
        Args:
            n_timesteps (int): Number of time steps to simulate
            params (dict): Dictionary of model parameters (A, Q, B, C, R, mu_0, sigma_0)
            options (dict): Options for LDS simulation
        Returns:
            tuple: Simulated state and observation arrays
        """
        n_dim_state = params['A'].shape[0]
        n_dim_obs = params['B'].shape[0]

        if options['inputs']:
            inp = options['inp']
        else:
            inp = np.zeros([n_timesteps, np.shape(params['C'])[1]])

        if 'h' not in params:
            params['h'] = np.zeros((1, n_dim_state))

        zi = scipy.stats.multivariate_normal(cov=params['Q'], allow_singular=True).rvs(n_timesteps)
        eta = scipy.stats.multivariate_normal(cov=params['R'], allow_singular=True).rvs(n_timesteps)

        state = np.zeros((n_timesteps, n_dim_state))
        obs = np.zeros((n_timesteps, n_dim_obs))

        for t in range(n_timesteps):
            if t == 0:
                state[t] = scipy.stats.multivariate_normal(mean=params['mu_0'] + params['C'] @ inp[0],
                                                     cov=params['sigma_0']).rvs(1)
            else:
                state[t] = params['A'] @ state[t - 1] + params['h'] + params['C'] @ inp[t] + zi[t]

            obs[t] = params['B'] @ state[t] + eta[t]

        return state, obs


    def eig_dynamics_features(self, matrix: np.array) -> tuple:
        eigvals, eigvecs = np.linalg.eig(matrix)
        idx_sorted = np.flip(np.argsort(eigvals))
        return eigvals[idx_sorted], eigvecs[:, idx_sorted]
        
    
    def svd_control_features(self, A: np.array, C: np.array) -> tuple:
        CC = control.ctrb(A, C)
        u,s,vh = np.linalg.svd(CC)
        return s,u

    def align_vectors(self, v1: np.array, v2: np.array) -> tuple:
        angle = np.dot(v1, v2)
        if angle < 0:
            v2 = -v2
        return v1, v2


    # Outlier Detection (IQR)
    def detect_ll_outlier(self, data: np.array) -> np.array:
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        ex = np.where(data <= (Q1-1.5*IQR))
        return ex


    def TwoSampleT2Test(self, X: np.array, Y: np.array) -> tuple:
            nx, p = X.shape
            ny, _ = Y.shape
            delta = np.mean(X, axis=0) - np.mean(Y, axis=0)
            Sx = np.cov(X, rowvar=False)
            Sy = np.cov(Y, rowvar=False)
            S_pooled = ((nx-1)*Sx + (ny-1)*Sy)/(nx+ny-2)
            t_squared = (nx*ny)/(nx+ny) * np.matmul(np.matmul(delta.transpose(), np.linalg.inv(S_pooled)), delta)
            statistic = t_squared * (nx+ny-p-1)/(p*(nx+ny-2))
            F = scipy.stats.f(p, nx+ny-p-1)
            p_value = 1 - F.cdf(statistic)
            return t_squared, statistic, p_value


    def group_difference_table(self, data: np.array, df: pd.DataFrame, stats_test: Optional[str]='TT', column_names: Optional[list]=None, index_names: Optional[list]=None, print_ms: Optional[bool]=True):
        idx = ~np.isnan(data).any(axis=0)
    
        group_data = {'group1': [], 'group2': []}
        means, stds = [], []
        for g in range(2):
            if data.ndim == 1:
                group_data['group' + str(g+1)] = data[(df['group']==g)&idx]
            else: 
                group_data['group' + str(g+1)] = data[:,(df['group']==g)&idx]
            means.append(np.round(np.nanmean(group_data['group' + str(g+1)].T,axis=0),2))
            stds.append(np.round(np.nanstd(group_data['group' + str(g+1)].T,axis=0),2))
    
        if stats_test == 'TT':
            [statistic,p] = scipy.stats.ttest_ind(group_data['group1'].T, \
                                          group_data['group2'].T)
    
        if stats_test == 'MWU':
                [statistic,p] = scipy.stats.mannwhitneyu(group_data['group1'].T, \
                                              group_data['group2'].T)
    
        if data.ndim == 1:
            df_stats = pd.concat((pd.DataFrame([str(means[j]) + ' ± ' + str(stds[j]) for j in range(2)]).T, \
                                                        pd.DataFrame([statistic, p]).T), axis=1)
        else: 
            df_stats = pd.concat((pd.DataFrame([[str(means[j][i]) + ' ± ' + str(stds[j][i]) for i in range(len(data))] for j in range(2)]).T, \
                                                        pd.DataFrame([statistic, p]).T), axis=1)
    
        if column_names != None: 
            df_stats.columns = column_names
            
        if index_names != None: 
            df_stats.index = index_names
            
        if print_ms: print('multiple comparison: p<' + str(0.05/data.ndim))
        
        if 'pvalue' in df_stats.columns: df_stats = self.replace_small_pvalues(df_stats)
        
        return df_stats


    def results_summary_to_dataframe(self, results) -> pd.DataFrame:
        '''take the result of an statsmodel results table and transforms it into a dataframe'''
        pvals = results.pvalues
        coeff = results.params
        conf_lower = results.conf_int()[0]
        conf_higher = results.conf_int()[1]
    
        results_df = pd.DataFrame({"pvals":pvals,
                                   "coeff":coeff,
                                   "conf_lower":conf_lower,
                                   "conf_higher":conf_higher
                                    })
    
        #Reordering...
        results_df = results_df[["coeff","conf_lower","conf_higher","pvals"]]
        return results_df


    def depression_relations(self, data: np.array, df: pd.DataFrame, column_names: Optional[list]=None, index_names: Optional[list]=None) -> pd.DataFrame:
    
        if data.ndim == 1: D = 1
        else: D, T = np.shape(data)
    
        tmp = pd.DataFrame()

        if D > 1:
            for i in range(D):
                idx = ~np.isnan(df) & ~np.isnan(data[i,:])
                [r,p] = scipy.stats.spearmanr(df[idx],data[i,idx])
                tmp = pd.concat((tmp,pd.DataFrame([np.round(r,2),p])),axis=1)
        else:
            idx = ~np.isnan(df) & ~np.isnan(data)
            [r,p] = scipy.stats.spearmanr(df[idx],data[idx])
            tmp = pd.DataFrame([np.round(r,2),p])

        tmp = tmp.T
        
        if column_names != None: 
            tmp.columns = column_names
            
        if index_names != None: 
            tmp.index = index_names
            
        if 'pvalue' in tmp.columns: tmp = self.replace_small_pvalues(tmp)
        
        return tmp


    def input_impact(self, df: pd.DataFrame) -> np.array:
        timing = df['timing'] - 1
        mood_data = df['data']
        inputs = df['inp']
    
        input_weights = np.full((4,inputs.shape[0]), np.nan)
        for i in range(inputs.shape[0]):
            input_idx = (inputs[i,timing]==1)
            input_timing = timing[input_idx]
            input_timing_pre = timing[np.append(input_idx[1:], False)].squeeze()
            if 0 in input_timing:
                change = mood_data[:,input_timing[1:]] - mood_data[:,input_timing_pre]
            else:
                change = mood_data[:,input_timing] - mood_data[:,input_timing_pre]
    
            input_weights[:,i] = np.mean(change, axis=1)
        return input_weights
    
        

