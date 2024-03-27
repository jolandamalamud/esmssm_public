import pandas as pd
import numpy as np
from numpy import ma
from statsmodels.tsa import stattools
from scipy.stats import skew,kurtosis
from typing import Optional
from ssm_stats import EsmssmStats
stats = EsmssmStats()

class EsmssmDataPrep:

    def __init__(self):
        pass
        
    
    def load_data(self, path: str, file: str, dataset: str) -> pd.DataFrame:
        
        if 'dta' in file:
            df = pd.read_stata(path + file, convert_categoricals=False)
        elif 'csv' in file:
            df = pd.read_csv(path + file, delimiter=';')
         
        df.name = dataset
        return df
          
    
    def time_to_minutes(self, df: pd.DataFrame) -> np.array:
        
        if df.name == 'leemput_patients':
            time0 = np.ceil(df['beeptime'][df['beeptime'] < 1] * 24 * 60)
            time = pd.concat([time0, df['beeptime'][df['beeptime'] > 1]])
            time = (np.floor(time / 60) * 60 + time % 60 + ((df['dayno'] - 1) * 24 * 60))

        elif df.name == 'leemput_controls':
            time = (np.floor(df['beeptime'] / 60) * 60 + df['beeptime'] % \
                    60 + ((df['dayno'] - 1) * 24 * 60))

        elif df.name == 'twindata':
            time = (np.floor(df['beeptime'] / 100) * 60 + df['beeptime'] % \
                    100 + ((df['dayno'] - 1) * 24 * 60))
            
        elif df.name == 'singledataset':
            time = pd.to_datetime(df['date'] + ' ' + df['beeptime'], format='%d/%m/%y %H:%M').astype(int)/(60*1e9)
            time = time - time[0]

        return time
    

    def psychiatric_score(self, df: pd.DataFrame, dataset: str) -> tuple:
        dep_score, labels = [], []
        if dataset == 'leemput_patients':
            # hdrs score at baseline - 6-8 week follow up hdrs score
            dep_score.append(df['hdrs1'].mean(skipna=True))
            dep_score.append(df['hdrs3'].mean(skipna=True))
            dep_score.append(dep_score[1] - dep_score[0])
            labels = ['dep1', 'dep2', 'change']
            
        elif dataset == 'leemput_controls': 
            # scl depression item score at baseline - follow up after 3 months
            dep_score.append(df['dep1'].mean(skipna=True))
            dep_score.append(df['dep2'].mean(skipna=True))
            dep_score.append(dep_score[1] - dep_score[0])
            labels = ['dep1', 'dep2', 'change']
            
        elif dataset == 'twindata':
            for i in ['dep_9f', 'anx_9f', 'tot']:
                dep_score.append(df['_1_scl90_' + i].mean(skipna=True))
                dep_score.append(df['_2_scl90_' + i].mean(skipna=True))
                dep_score.append(dep_score[1] - dep_score[0])
                labels = labels + [i + '1', i + '2', i + '_change']
        
        return dep_score, labels
    
    
    def creating_df(self, df: pd.DataFrame, columns: list, labels: list, dataset: str, inputs: Optional[bool]=False, input_columns: Optional[bool]=None, input_labels: Optional[bool]=None) -> tuple:
        df.name = dataset
        df['originalid'] = df['subjno'].copy()
        df['subjno'] = df['subjno'].astype('category').cat.codes + 1
        df['timing'] = self.time_to_minutes(df)

        exclusion_code = []
        new_df = pd.DataFrame()
        df_inputs_all = pd.DataFrame()

        # subject info
        n_sub = len(np.unique(df['subjno'][~np.isnan(df['subjno'])]))
        sub_list = np.unique(df['subjno'][~np.isnan(df['subjno'])])

        print('prep data for ' + df.name + ' (N=' + str(n_sub) + '):', end='\n' + 100 * '-' + '\n')

        iid = 1
        for sub_id in sub_list:
            print('setting up mood time series for subject {}'.format(sub_id))
            idx_sub = df['subjno'] == sub_id

            df_subject = df.loc[idx_sub].reset_index(drop=True)
            T = len(df_subject)

            timing = df_subject['timing'].reset_index(drop=True)
            mood = df_subject[columns].reset_index(drop=True)
            nan_idx = mood.isna().any(axis=1) + timing.isna()

            # exclude if no data at all
            if nan_idx.all():
                print('deleting subject {}, because no data'.format(sub_id))
                exclusion_code.append([sub_id, 1])
                continue

            timing = timing[~nan_idx].to_numpy().astype(int)
            mood = mood[~nan_idx].to_numpy().astype(int).T

            # exclude if esm less than XX days and less than XX timepoints
            n_days = np.max(df_subject['dayno']);
            if n_days < 5:
                print('deleting subject {}, only {} days of esm'.format(sub_id, n_days))
                exclusion_code.append([sub_id, 2])
                continue

            if len(mood[0]) < 40:
                print('deleting subject {}, because compliance = {}'.format(sub_id, len(mood[0])))
                exclusion_code.append([sub_id, 3])
                continue    

            # exclude if no variability in mood items
            if all(np.diag(np.cov(mood)) == 0):
                print('deleting subject {}, because no variability in mood'.format(sub_id))
                exclusion_code.append([sub_id, 4])
                continue
            
            # check time points
            if not (np.sort(timing) == timing).all():
                print('deleting subject {}, because something wrong with time'.format(sub_id))
                exclusion_code.append([sub_id, 5])
                continue
             

            df_subject_final = pd.DataFrame(mood.T, \
                                            columns=['mood_' + i for i in labels])
            df_subject_final['timing'] = timing - timing[0] + 1
            df_subject_final.insert(0,'id',iid)
            df_subject_final.insert(1,'subid',sub_id)
            df_subject_final.insert(2,'originalid', df_subject['originalid'].iloc[0])
            
            score, score_labels = self.psychiatric_score(df_subject, df.name)
            for i, l in enumerate(score_labels):
                df_subject_final[l] = score[i]

            new_df = pd.concat((new_df, df_subject_final),axis=0)
            
            # create input matrix
            if inputs:
                df_input = self.create_input_matrix(df_subject[~nan_idx], input_columns, input_labels)
                inp_nocoll,input_idx = self.remove_collinear_rows(df_input.to_numpy())
                df_input_nocoll = pd.DataFrame(inp_nocoll, columns=df_input.columns[input_idx])
                df_input_nocoll.insert(0,'id',iid)
                df_input_nocoll.insert(1, 'originalid', sub_id)
                df_input_nocoll[df_input_nocoll>1] = 1

                df_inputs_all = pd.concat((df_inputs_all, df_input_nocoll))
            
            iid += 1 

        return new_df, pd.DataFrame(exclusion_code, columns=['subid', 'exclusion_code']), df_inputs_all

    
    def remove_collinear_rows(self, inp: np.array) -> tuple:
        tolerance = 1e-10
        q, r = np.linalg.qr(inp)
        inp_fullrank = inp[:, np.abs(np.diag(r)) > tolerance]
        idx_inputs_kept = np.arange(np.shape(inp)[1])[np.abs(np.diag(r)) > tolerance]
        
        return inp_fullrank, idx_inputs_kept

    
    def create_input_matrix(self, df: pd.DataFrame, columns: list, labels: list) -> pd.DataFrame:
        input_df = pd.DataFrame()
        for k,i in enumerate(columns):
            specific_df = pd.DataFrame()
            for j in i:
                subspecific_df = pd.get_dummies(df[j])
                for c in subspecific_df.columns:
                    if c in labels[k].keys():
                        subspecific_df = subspecific_df.rename(columns={c: labels[k][c]})
                    else:
                        subspecific_df = subspecific_df.drop(columns=c)
                specific_df = specific_df.add(subspecific_df, fill_value=0)
            input_df = pd.concat((input_df, specific_df), axis=1)

        return input_df.fillna(0)
    
    
    def creating_input_df(self, df: pd.DataFrame, sub_list: list, columns: list, labels: list) -> pd.DataFrame:
        df_inputs_all = pd.DataFrame()
        iid = 1
        for sub_id in sub_list: 
        #     print('setting up input time series for subject {}'.format(sub_id))
            idx_sub = df['subjno'] == sub_id
            df_subject = df.loc[idx_sub].reset_index(drop=True)

            # create input matrix
            df_input = self.create_input_matrix(df_subject, columns, labels)
            inp_nocoll,input_idx = self.remove_collinear_rows(df_input.to_numpy())
            df_input_nocoll = pd.DataFrame(inp_nocoll, columns=df_input.columns[input_idx])
            df_input_nocoll.insert(0,'id',iid)
            df_input_nocoll.insert(1, 'originalid', sub_id)
            df_input_nocoll[df_input_nocoll>1] = 1

            df_inputs_all = pd.concat((df_inputs_all, df_input_nocoll))
            iid += 1
            
        return df_inputs_all
    

    def bimodality_coef(self, data: np.array) -> float:
        N = len(data)
        NN = (N-1)**2 / ((N-2)*(N-3))
        if np.nanvar(data) == 0: return np.nan
        s = skew(data, axis=0, bias=True)
        k = kurtosis(data, axis=0, bias=True)
        return (s**2 + 1) / (k + 3 * NN)
    
    
    def basic_mood_dynamics(self, data: np.array) -> dict:
        ns, T = data.shape
        basics = dict()
        basics['mean'] = np.nanmean(data,axis=1)
        basics['variance'] = np.nanvar(data,axis=1)
        basics['covariance'] = np.cov(data[:,~np.isnan(data[0,:])])
        basics['RMSSD'] = 1/(T - 1) * np.nansum(data[:,1:]-data[:,:-1],axis=1)
        autocorr_coeff, bimod_coeff = [], []
        for i in range(ns):
            if basics['variance'][i] != 0:
                autocorr = stattools.acf(data[i,:])
                autocorr_coeff.append(autocorr[1])
            else: autocorr_coeff.append(np.nan)
            bimod_coeff.append(self.bimodality_coef(data[i,:]))
        basics['AR'] = autocorr_coeff
        basics['bc'] = bimod_coeff
        return basics

    
    def mergeDictionary(self, dict_1: dict, dict_2: dict) -> dict:
        dict_3 = {**dict_1, **dict_2}
        for key, value in dict_3.items():
            if key in dict_1 and key in dict_2:
                dict_3[key] = np.dstack((dict_1[key] , dict_2[key]))
        return dict_3


    def extract_dynamic_features(self, data: np.array, control_features: Optional[bool]=False, input_columns: Optional[list]=None, inputs: Optional[np.array]=None):
        Nsj = len(data)
        m,n = np.shape(data[0]['A'])
        dynamics = dict()
        dynamics['matrix'], dynamics['matrix_hour'], dynamics['vec'], dynamics['vec_imag'], dynamics['VAR'], \
                dynamics['vec_aligned'], dynamics['convec'], dynamics['convec_aligned'] = [np.full((m,n,Nsj), np.nan) for _ in range(8)]
        dynamics['bias'], dynamics['val'], dynamics['val_imag'], \
            dynamics['diff'], dynamics['FP'], dynamics['conval']  = [np.full((m,Nsj), np.nan) for _ in range(6)]
        dynamics['stab'] = np.full((Nsj), np.nan)

        if control_features:
            ni = len(input_columns)
            dynamics['conmatrix'] = np.full((n,ni,Nsj), np.nan)
    
        basics = []
    
        for i in range(Nsj):
            dynamics['matrix'][:,:,i] = data[i]['A'] # minute matrix
            dynamics['stab'][i] = np.linalg.det(data[i]['A'])
            dynamics['matrix_hour'][:,:,i] = \
                np.linalg.matrix_power(data[i]['A'],60) # hourly matrix
            dynamics['bias'][:,i] = np.squeeze(data[i]['bias'])
            [e,v] = stats.eig_dynamics_features(data[i]['A'])
            dynamics['val'][:,i] = np.real(e)
            dynamics['vec'][:,:, i] = np.real(v)
            dynamics['diff'][:,i] = np.abs(np.nanmean(np.real(v)[:2,:],axis=0) \
                                           - np.nanmean(np.real(v)[2:,:],axis=0))
            if control_features:
                dynamics['conmatrix'][:,inputs[:,i]==1,i] = data[i]['C']
                [s,u] = stats.svd_control_features(data[i]['A'], data[i]['C'])
                dynamics['conval'][:,i] = np.real(s)
                dynamics['convec'][:,:, i] = np.real(u)
                
            basics.append(self.basic_mood_dynamics(np.squeeze(data[i]['data'][:,data[i]['timing']-1])))
            
            for k in range(m):
                foo, dynamics['vec_aligned'][:,k, i] = \
                stats.align_vectors([1,1,-1,-1],np.real(v[:,k]))
                if control_features:
                    foo, dynamics['convec_aligned'][:,k, i] = \
                    stats.align_vectors([1,1,-1,-1],np.real(u[:,k]))
        
            try:
                varresults, foo = stats.fit_simulate_var(np.squeeze(data[i]['data'][:,data[i]['timing']-1]).T)
                dynamics['VAR'][:,:,i] = varresults.coefs 
            except:
                pass
                
        basic_char = dict()
        for i in range(Nsj):
            basic_char = self.mergeDictionary(basic_char, basics[i])
        for key in basic_char.keys():
            basic_char[key] = np.squeeze(basic_char[key]) 

        return dynamics, basic_char
            

   
