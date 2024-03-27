import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import SubplotSpec
from typing import Optional

class EsmssmPlotting:

    def __init__(self):
        pass
    
    
    def pval_to_asterix(self, p: float, multiple_comparison: Optional[float]=None) -> str:
        if multiple_comparison == None: multiple_comparison = 1
        if p <= 0.001/multiple_comparison: asterix = '***'
        elif p <= 0.01/multiple_comparison: asterix = '**'
        elif p <= 0.05/multiple_comparison: asterix = '*'
        else: asterix = '$^{ns}$'
        return asterix

    
    def create_subtitle(self, fig: plt.Figure, grid: SubplotSpec, title: str, x_loc: float, y_loc: float):
        "Sign sets of subplots with title"
        # Use loc parameter to control title location
        row = fig.add_subplot(grid)
        # the '\n' is important
        row.set_title(f'{title}\n', fontweight='semibold', fontsize=28, x=x_loc, y=y_loc)
        # hide subplot
        row.set_frame_on(False)
        row.axis('off')

            
    def create_annotation(self, data: np.array, pval: Optional[float]=None, parametric: Optional[bool]=None) -> np.array:
        dim = np.shape(data)
        data = np.real(data)
        if len(dim) == 3:
            annotation = np.full([dim[0],dim[1]], np.nan, dtype=object)
            for i in range(dim[0]):
                for j in range(dim[1]):
                    if parametric is not None and parametric is True:
                        annotation[i,j] = str(np.round(np.nanmean(data[i,j,:]),3)) + '\n± ' \
                                        + str(np.round(np.nanstd(data[i,j,:])/np.sqrt(len(data[i,j,:])),2))
                    else:
                        annotation[i,j] = str(np.round(np.nanmedian(data[i,j,:]),3)) + '\n± ' \
                                        + str(np.round(scipy.stats.iqr(data[i,j,:], nan_policy='omit'),2))
                            
                    if pval is not None:
                        if pval[i,j] < 0.05:
                            annotation[i,j] = self.pval_to_asterix(pval[i,j]) + '\n' + annotation[i,j]

        elif len(dim) == 2:     
            annotation = np.full(dim[0], np.nan, dtype=object)
            for i in range(dim[0]):
                if parametric is not None and parametric is True:
                    annotation[i,j] = str(np.round(np.nanmean(data[i,j,:]),3)) + '\n± ' \
                                    + str(np.round(np.nanstd(data[i,j,:]),2))
                else:
                    annotation[i,j] = str(np.round(np.nanmedian(data[i,j,:]),3)) + '\n± ' \
                                    + str(np.round(scipy.stats.iqr(data[i,j,:], nan_policy='omit'),2))
                if pval is not None:
                    if pval[i] < 0.05:
                        annotation[i] = self.pval_to_asterix(pval[i]) + '\n' + annotation[i]

        return annotation  
    
    
    def plot_dynamics_matrices_group(self, matrices: np.array, group: list, pval: list, opt: dict):
        if 'fig' in opt:
            fig = opt['fig']
            ax = opt['ax']
        else:
            fig,ax = plt.subplots(1,2, figsize=(20,6))

        data_plot, annotation, vmin, vmax = [[] for i in range(4)]
        for i in range(2):
            data_plot.append(matrices[:,:,group==i])
            annotation.append(self.create_annotation(data_plot[i], pval.reshape(4,4), opt['parameteric']))
            if opt['parameteric'] is not None and opt['parameteric'] is True:
                vmin.append(np.min(np.nanmean(data_plot[i],axis=2)))
                vmax.append(np.max(np.nanmean(data_plot[i],axis=2)))
            else:
                vmin.append(np.min(np.nanmedian(data_plot[i],axis=2)))
                vmax.append(np.max(np.nanmedian(data_plot[i],axis=2)))
        for i in range(2):
            if opt['parameteric'] is not None and opt['parameteric'] is True: 
                d = np.nanmean(data_plot[i],axis=2)
            else: 
                d = np.nanmedian(data_plot[i],axis=2)
            ss = sns.heatmap(d, ax=ax[i], annot=annotation[i], fmt='', \
                             annot_kws={'color':'black'}, vmin = np.min(vmin), vmax=np.max(vmax), \
                             yticklabels=opt['xlabels'], xticklabels=opt['xlabels'], cmap=opt['cmap'][i], \
                             cbar=False, center=0)
            ss.set_yticklabels(labels=ss.get_yticklabels(), va='center')
            ax[i].set_title(opt['group_label'][i])
            sig = np.where(pval.reshape(4,4)*16 <= (0.05))
            for j in range(len(sig[0])):
                ss.add_patch(Rectangle((sig[1][j], sig[0][j]), 1, 1, ec='black', fc='none', lw=5))

                
    def boxplot_group(self, data: np.array, group: list, pval: list, opt: dict, annot: Optional[bool]=None, msc: Optional[bool]=None):
        if 'fig' in opt:
            fig = opt['fig']
            ax = opt['ax']
        else:
            fig, axs = plt.subplots(1,2, figsize=(20,6))

        df_melt = pd.DataFrame(group)
        df_melt = pd.concat((df_melt,pd.DataFrame(data, columns=opt['xticklabels'])),axis=1)
        df_melt = pd.melt(df_melt, id_vars='group', value_vars=opt['xticklabels'], var_name='mood', value_name='value')

        hue_plot_params = {
            'data': df_melt,
            'x': 'mood',
            'y': 'value',
            "hue": "group",
            "palette": opt['color'],
            "showmeans": True,    }

        b = sns.boxplot(**hue_plot_params, ax=ax, meanprops={"markerfacecolor":"white", "markeredgecolor":"black"})
        b.set(xlabel='', ylabel=opt['ylabel'])

        pairs = []
        for i in opt['xticklabels']:
            pairs.append([(i, 0), (i, 1)])

        # Add annotations
        if annot:
            annotator = Annotator(b, pairs, **hue_plot_params);
            annotator.configure(verbose=0, test='t-test_welch', comparisons_correction=msc)
            annotator.apply_and_annotate()
        return b


    def covariance_plot(self, c: np.array, fig: plt.Figure, outer: SubplotSpec, labels: list):
        n,m,Nsj = np.shape(c['true'])
        for i in range(m):
            for j in range(m):
                ax = plt.Subplot(fig, outer[i,j])
                ax = fig.add_subplot(ax)
                sns.regplot(x=c['true'][i,j,:], y=c['var'][i,j,:],ax=ax, truncate=False, scatter_kws={'s':5, 'marker':'+'}, color='gold')
                sns.regplot(x=c['true'][i,j,:], y=c['kf'][i,j,:],ax=ax, truncate=False, scatter_kws={'s':5, 'marker':'o'})
                ax.set(xticks=[], yticks=[])    
                if j == 0: ax.set(ylabel=labels[i])
                if i == m-1: ax.set(xlabel=labels[j])


    def ar_plot(self, ar_lags: dict, fig: plt.Figure, outer: SubplotSpec, labels: list):
        axs = None
        for i in range(len(labels)):
            ax = plt.Subplot(fig, outer[-1,i])
            ax = fig.add_subplot(ax, sharey=axs, sharex=axs)
            for j, k in enumerate(ar_lags.keys()):
                y = np.nanmean((np.array(ar_lags[k])[:,i,:10]),axis=0)
                x = range(len(y))
                error = np.nanstd((np.array(ar_lags[k])[:,i,:10]),axis=0)/np.sqrt(np.shape(ar_lags[k])[0])
                ax.plot(x, y, ['r-o', 'b-o', 'y-o'][j], label=['observations', 'KF stimulation', 'VAR simulation'][j])
                ax.fill_between(x, y-error, y+error, alpha=0.2, color=['red', 'blue', 'yellow'][j])
            ax.set(xlabel='lags')
            ax.set_title(labels[i])
            if i == 0: ax.set(ylabel='autocorrelation coefficient')
            if i != 0: ax.set(yticklabels=[])
            if i == 3: ax.legend()
            axs = ax
            

    def color_shading(self, colors: list):
        c = [colors[0],colors[0], "white", colors[1], colors[1]]
        v = [0,0.25,.5,0.75,1]
        l = list(zip(v,c))
        return LinearSegmentedColormap.from_list('rg',l, N=256)

    
    def color_definition(self) -> tuple:
        # colormap for colorblind
        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']
        k = sns.color_palette(CB_color_cycle)
        group_color = [CB_color_cycle[0], CB_color_cycle[1]]
        negpos_color = [CB_color_cycle[7], CB_color_cycle[2]]
        minusplus_color = [CB_color_cycle[6], CB_color_cycle[5]]
        patients_shading = self.color_shading([CB_color_cycle[6], CB_color_cycle[0]])
        controls_shading = self.color_shading([CB_color_cycle[6], CB_color_cycle[1]])

        return group_color, negpos_color, minusplus_color, patients_shading, controls_shading

