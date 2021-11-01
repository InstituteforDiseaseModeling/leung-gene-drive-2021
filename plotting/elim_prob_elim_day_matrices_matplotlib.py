from matplotlib import rcParams
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

params = {'axes.labelsize': 10,
          'axes.titlesize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10}
rcParams.update(params)

rcParams['pdf.fonttype'] = 42

##
# -------- Set figures to plot
plot_elim_probs = 1
anno_ep_all_cols = 1  # 0 = don't annotate, 1 = annotate bottom square

plot_elim_days = 1

outer_axes = 'rc_dord1'  # choose: 'rc_dord1', 'rr0orrr20_sneorse2'

##
# -------- Define fxns and colormaps

# - Define colors maps
# --> Plotly colors
# import plotly.colors as colors
# greens_full = colors.get_colorscale('greens')
# greens = greens_full[1:]
# for i in range(0, len(greens)):
#     greens[i][0] = i / (len(greens) - 1)
# --> The following uses the rgb values from "greens" above
greensnow = [(229, 245, 224), (199, 233, 192),
             (161, 217, 155), (116, 196, 118),
             (65, 171, 93), (35, 139, 69),
             (0, 109, 44), (0, 68, 27)]
greensnow = [tuple(i/255 for i in tuplenow) for tuplenow in greensnow]
greens = LinearSegmentedColormap.from_list(name='cmapnow', colors=greensnow, N=256)
greens_rnow = greensnow[::-1]
greens_r = LinearSegmentedColormap.from_list(name='cmaprnow', colors=greensnow, N=256)

# - Define characteristics for matrix text annotations
anno_text_colors = ["black", "white"]

# - Define characteristics and fxns for symbol annotations
anno_symbol_size = 10
anno_symbol_color = 'orange'
anno_y_offset = -0.67


def monotonic_increasing(x):
    if len(np.unique(x)) > 1:
        dx = np.diff(x)
        return np.all(dx >= 0)
    else:
        return False


def monotonic_decreasing(x):
    if len(np.unique(x)) > 1:
        dx = np.diff(x)
        return np.all(dx <= 0)
    else:
        return False


def wiggly(x):
    dx = np.diff(x)
    if (len(np.unique(x)) > 1) & \
            ~(np.all(dx <= 0)) & \
            ~(np.all(dx >= 0)):
        return True

##
# -------- Setup datasets and sweep var values
wi_names_ls = [
    'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne',
    'spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2',
    'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne',
    'spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2',
    'spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne',
    'spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2',
    'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne',
    'spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2',
    'spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne',
    'spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2'
]
num_sweep_vars_ls = [
    4,
    4,
    4,
    4,
    4, 4, 4, 4, 4, 4
]
drive_types_ls = [
    'classic',
    'integral',
    'classic',
    'integral',
    'classic', 'integral', 'classic', 'integral', 'classic', 'integral'
]
data_dir = '..\\csvs'
fig_dir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\presentations_writeups\\' \
          'gene_drive_paper\\figures\\elim_prob_day_matrices\\' + outer_axes + '_outside'
num_seeds = 20  # num of seeds per sim

for iwi, wi_name in enumerate(wi_names_ls):
    drive_type = drive_types_ls[iwi]
    num_sweep_vars = num_sweep_vars_ls[iwi]

    if num_sweep_vars == 6:
        if drive_type == 'classic':
            allvardefs = {'rc': 1, 'd': 1,
                          'rr0': 0, 'sne': 0,
                          'rd': 180, 'nn': 6}
            allvarvals = {
                'rc': [0.5, 0.6, 0.7, 0.8, 0.9, 1],
                'd': [0.9, 0.95, 1],
                'rr0': [0, 0.001, 0.01, 0.1],
                'sne': [0, 0.1, 0.2, 0.3, 0.4, 0.5],
                'rd': [180, 240, 300, 360, 420, 480, 545],
                'nn': [6, 12]
            }
    elif num_sweep_vars == 4:
        if drive_type == 'classic':
            allvardefs = {'rc': 1, 'd': 1,
                          'sne': 0, 'rr0': 0}
            allvarvals = {
                'rc': [0.5, 0.6, 0.7, 0.8, 0.9, 1],
                'd': [0.9, 0.95, 1],
                'rr0': [0, 0.001, 0.01, 0.1],
                'sne': [0, 0.1, 0.2, 0.3, 0.4, 0.5]
            }
        elif drive_type == 'integral':
            allvardefs = {'rc': 1, 'd1': 1,
                          'se2': 0, 'rr20': 0}
            allvarvals = {
                'rc': [0.5, 0.6, 0.7, 0.8, 0.9, 1],
                'd1': [0.9, 0.95, 1],
                'rr20': [0, 0.001, 0.01, 0.1],
                'se2': [0, 0.1, 0.2, 0.3, 0.4, 0.5]
            }

    ##
    # -------- Load data
    dfe = pd.read_csv(os.path.join(data_dir, 'dfe_' + wi_name + '.csv'))
    dfed = pd.read_csv(os.path.join(data_dir, 'dfed_' + wi_name + '.csv'))

    if num_sweep_vars == 6:
        dfe['rd'] = dfe['rd'].fillna(allvardefs['rd'])
        dfe['nn'] = dfe['nn'].fillna(allvardefs['nn'])
        dfe = dfe[dfe['rd'] == allvardefs['rd']]
        dfe = dfe[dfe['nn'] == allvardefs['nn']]
        dfed['rd'] = dfed['rd'].fillna(allvardefs['rd'])
        dfed['nn'] = dfed['nn'].fillna(allvardefs['nn'])
        dfed = dfed[dfed['rd'] == allvardefs['rd']]
        dfed = dfed[dfed['nn'] == allvardefs['nn']]

    ##
    # -------- Set up x/y axes
    # - Set matrix x/y vars, overall x/y vars
    if drive_type == 'classic':
        if outer_axes == 'rc_dord1':
            # - rc/d on outer axes
            mat_xvar = 'rr0'
            mat_yvar = 'sne'
            ov_xvar = 'rc'
            ov_yvar = 'd'
        elif outer_axes == 'rr0orrr20_sneorse2':
            # - rr0/sne on outer axes
            mat_xvar = 'rc'
            mat_yvar = 'd'
            ov_xvar = 'sne'
            ov_yvar = 'rr0'
    elif drive_type == 'integral':
        if outer_axes == 'rc_dord1':
            # - rc/d1 on outer axes
            mat_xvar = 'rr20'
            mat_yvar = 'se2'
            ov_xvar = 'rc'
            ov_yvar = 'd1'
        elif outer_axes == 'rr0orrr20_sneorse2':
            # - rr0/sne on outer axes
            mat_xvar = 'rc'
            mat_yvar = 'd1'
            ov_xvar = 'se2'
            ov_yvar = 'rr20'

    # - Set variable title strings
    ov_xvar_strnow = ov_xvar + ' = '
    if ov_xvar == 'rc':
        ov_xvar_titlestrnow = 'Transmission-blocking effectiveness (' + ov_xvar + ')'
    elif (ov_xvar == 'd') or (ov_xvar == 'd1'):
        ov_xvar_titlestrnow = 'Drive efficiency (' + ov_xvar + ')'
    elif (ov_xvar == 'sne') or (ov_xvar == 'se2'):
        ov_xvar_titlestrnow = 'Fitness cost (' + ov_xvar + ')'
    elif (ov_xvar == 'rr0') or (ov_xvar == 'rr20'):
        ov_xvar_strnow = 'Pre-existing resistance (' + ov_xvar + ') = '

    if ov_yvar == 'rc':
        ov_yvar_strnow = 'Transm.-blocking\neffectiv. (' + ov_yvar + ') = '
    elif (ov_yvar == 'd') or (ov_yvar == 'd1'):
        ov_yvar_strnow = 'Drive efficiency\n(' + ov_yvar + ') = '
    elif (ov_yvar == 'sne') or (ov_yvar == 'se2'):
        ov_yvar_strnow = 'Fitness cost\n(' + ov_yvar + ') = '
    elif (ov_yvar == 'rr0') or (ov_yvar == 'rr20'):
        ov_yvar_strnow = 'Initial resistance\n (' + ov_yvar + ') = '
    else:
        ov_yvar_strnow = ov_yvar + ' ='

    if mat_xvar == 'rc':
        mat_xvar_strnow = 'Transmission-blocking effectiv. (' + mat_xvar + ')'
    elif (mat_xvar == 'd') or (mat_xvar == 'd1'):
        mat_xvar_strnow = 'Drive efficiency (' + mat_xvar + ')'
    elif (mat_xvar == 'sne') or (mat_xvar == 'se2'):
        mat_xvar_strnow = 'Fitness cost (' + mat_xvar + ')'
    elif (mat_xvar == 'rr0') or (mat_xvar == 'rr20'):
        mat_xvar_strnow = 'Initial resistance (' + mat_xvar + ')'
    else:
        mat_xvar_strnow = mat_xvar

    if mat_yvar == 'rc':
        mat_yvar_strnow = 'Transmission-blocking effectiv. (' + mat_yvar + ')'
    elif (mat_yvar == 'd') or (mat_yvar == 'd1'):
        mat_yvar_strnow = 'Drive efficiency (' + mat_yvar + ')'
    elif (mat_yvar == 'sne') or (mat_yvar == 'se2'):
        mat_yvar_strnow = 'Fitness cost (' + mat_yvar + ')'
    elif (mat_yvar == 'rr0') or (mat_yvar == 'rr20'):
        mat_yvar_strnow = 'Initial resistance (' + mat_yvar + ')'
    else:
        mat_yvar_strnow = mat_yvar

    # - Set over variable values
    ov_xvar_vals = allvarvals[ov_xvar]
    ov_yvar_vals = allvarvals[ov_yvar]

    ##
    # -------- Create elim prob matrix
    if plot_elim_probs == 1:

        dfesm = dfe[dfe[mat_xvar].isin(allvarvals[mat_xvar]) &
                    dfe[mat_yvar].isin(allvarvals[mat_yvar])]

        fig, axes = plt.subplots(nrows=len(ov_yvar_vals), ncols=len(ov_xvar_vals),
                                 figsize=(11, 6), sharex=True, sharey=True)
        fig.subplots_adjust(hspace=0.07, wspace=0.09)  # adjust spacing btwn subplots

        for iov_yvar, ov_yvar_val in enumerate(ov_yvar_vals):
            for iov_xvar, ov_xvar_val in enumerate(ov_xvar_vals):

                # - Compute heatmap values
                allvardefsnow = {k: v for k, v in allvardefs.items() if k not in [mat_xvar, mat_yvar, ov_xvar, ov_yvar]}
                dfenow = dfesm
                if len(allvardefsnow) > 0:
                    for k, v in allvardefsnow.items():
                        dfenow = dfenow[dfenow[k] == v]
                        dfenow.drop(columns=[k], inplace=True)
                dfenow = dfenow[dfenow[ov_xvar] == ov_xvar_val]
                dfenow = dfenow[dfenow[ov_yvar] == ov_yvar_val]
                dfenow.drop(columns=[ov_xvar, ov_yvar], inplace=True)
                dfenownow = (dfenow.groupby([mat_xvar, mat_yvar])[
                                 'True_Prevalence_elim'].sum() / num_seeds).reset_index()
                matnow = dfenownow.pivot_table(index=[mat_yvar], columns=[mat_xvar], values='True_Prevalence_elim')
                matnow = np.flipud(matnow.values)

                # - Plot heat map
                ax = axes[iov_yvar][iov_xvar]
                im = ax.imshow(matnow, cmap=greens, vmin=0, vmax=1, aspect="auto")  # aspect auto needed for wspace
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['left'].set_visible(False)
                ax.tick_params(axis='both', length=2, which='major')
                ax.tick_params(axis='x', labelrotation=-45)

                # - Set ticks, ticklabels, overall x and y var labels
                if iov_yvar == 0:
                    # Column titles (overall x labels)
                    ax.set_title(ov_xvar_strnow + str(ov_xvar_val))
                    ax.axes.xaxis.set_visible(False)
                elif iov_yvar == len(ov_yvar_vals)-1:
                    ax.set_xticks(np.arange(len(allvarvals[mat_xvar])))
                    ax.set_xticklabels(allvarvals[mat_xvar])
                else:
                    ax.axes.xaxis.set_visible(False)

                if iov_xvar == 0:
                    ax.set_yticks(np.arange(len(allvarvals[mat_yvar])))
                    ax.set_yticklabels(allvarvals[mat_yvar][::-1])
                    # Row titles (overall y labels)
                    ax.set_ylabel(ov_yvar_strnow + str(ov_yvar_val), rotation=-90, labelpad=542)
                    ax.yaxis.set_label_position("right")
                else:
                    ax.axes.yaxis.set_visible(False)

                # - Set matrix x and y var labels
                fig.text(0.5, 0.015, mat_xvar_strnow, ha='center', va='center')
                fig.text(0.085, 0.5, mat_yvar_strnow, ha='center', va='center', rotation='vertical')
                fig.text(0.085, 0.5, mat_yvar_strnow, ha='center', va='center', rotation='vertical')

                # - Add elim prob annotations
                anno_textcol_threshold = 0.5
                for i in range(len(allvarvals[mat_yvar])):
                    for j in range(len(allvarvals[mat_xvar])):
                        colornow = anno_text_colors[int(im.norm(matnow[i, j]) >= anno_textcol_threshold)]
                        annonow = round(round(matnow[i, j]/0.05)*0.05, 2)
                        text = ax.text(j, i, str(annonow), ha="center", va="center",
                                       color=colornow, fontsize=8)

                # - Categorize columns into monotonic increasing, mt decreasing, wiggly
                incr_cols = []
                decr_cols = []
                wiggly_cols = []
                for icol in range(0, matnow.shape[1]):
                    colnow = np.flip(matnow[:, icol])
                    incr_cols.append(monotonic_increasing(colnow))
                    decr_cols.append(monotonic_decreasing(colnow))
                    wiggly_cols.append(wiggly(colnow))

                # - Annotate columns (mt increasing, mt decreasing, wiggly)
                mono_decr_xs = np.where(decr_cols)[0]
                mono_incr_xs = np.where(incr_cols)[0]
                wiggly_xs = np.where(wiggly_cols)[0]
                if anno_ep_all_cols == 1:
                    if len(mono_decr_xs) > 0:
                        for anno_x in mono_decr_xs:
                            text = ax.text(anno_x, len(allvarvals[mat_yvar]) + anno_y_offset,
                                           r'$\blacktriangledown$', ha="center", va="center",
                                           color=anno_symbol_color, fontsize=anno_symbol_size)
                    if len(mono_incr_xs) > 0:
                        for anno_x in mono_incr_xs:
                            text = ax.text(anno_x, len(allvarvals[mat_yvar]) + anno_y_offset,
                                           r'$\blacktriangle$', ha="center", va="center",
                                           color=anno_symbol_color, fontsize=anno_symbol_size)
                    if len(wiggly_xs) > 0:
                        for anno_x in wiggly_xs:
                            text = ax.text(anno_x, len(allvarvals[mat_yvar]) + anno_y_offset,
                                           '~', ha="center", va="center",
                                           color=anno_symbol_color, fontsize=anno_symbol_size, fontweight='bold')

        # - Add colorbar
        fig.subplots_adjust(right=0.9)
        cbar_ax = fig.add_axes([0.95, 0.1, 0.02, 0.8])  # left, bottom, width, height
        fig.colorbar(im, cax=cbar_ax)

        # - Save fig
        fig_file_png = os.path.join(fig_dir, wi_name + '_elim_probs.png')
        fig_file_pdf = os.path.join(fig_dir, wi_name + '_elim_probs.pdf')
        plt.savefig(fig_file_png, bbox_inches="tight", dpi=300)
        plt.savefig(fig_file_pdf, bbox_inches="tight", dpi=300)
        plt.show()

    ##
    # -------- Create elim day matrix
    if plot_elim_days == 1:

        dfedsm = dfed[dfed[mat_xvar].isin(allvarvals[mat_xvar]) &
                      dfed[mat_yvar].isin(allvarvals[mat_yvar])]

        fig, axes = plt.subplots(nrows=len(ov_yvar_vals), ncols=len(ov_xvar_vals),
                                 figsize=(11, 6), sharex=True, sharey=True)
        fig.subplots_adjust(hspace=0.07, wspace=0.09)  # adjust spacing btwn subplots

        for iov_yvar, ov_yvar_val in enumerate(ov_yvar_vals):
            for iov_xvar, ov_xvar_val in enumerate(ov_xvar_vals):

                # - Compute heatmap values
                allvardefsnow = {k: v for k, v in allvardefs.items() if k not in [mat_xvar, mat_yvar, ov_xvar, ov_yvar]}
                dfednow = dfedsm
                if len(allvardefsnow) > 0:
                    for k, v in allvardefsnow.items():
                        dfednow = dfednow[dfednow[k] == v]
                        dfednow.drop(columns=[k], inplace=True)
                dfednow = dfednow[dfednow[ov_xvar] == ov_xvar_val]
                dfednow = dfednow[dfednow[ov_yvar] == ov_yvar_val]
                dfednow.drop(columns=[ov_xvar, ov_yvar], inplace=True)
                dfednow.loc[dfednow['True_Prevalence_elim'] == False,
                            'True_Prevalence_elim_day'] = np.nan
                dfednow.drop(columns=['True_Prevalence_elim'], inplace=True)
                dfednownow = (dfednow.groupby([mat_xvar, mat_yvar])['True_Prevalence_elim_day'].mean()).reset_index()
                matnow = dfednownow.pivot_table(index=[mat_yvar], columns=[mat_xvar],
                                                values='True_Prevalence_elim_day', dropna=False)
                matnow = (matnow / 365)
                matnow = np.flipud(matnow.values)

                # - Plot heat map
                ax = axes[iov_yvar][iov_xvar]
                greens.set_bad('gainsboro')
                im = ax.imshow(matnow, cmap=greens, vmin=(dfed['True_Prevalence_elim_day'] / 365).min(),
                               vmax=(dfed['True_Prevalence_elim_day'] / 365).max(), aspect="auto")  # aspect auto needed for wspace
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['left'].set_visible(False)
                ax.tick_params(axis='both', length=2, which='major')
                ax.tick_params(axis='x', labelrotation=-45)

                # - Set ticks, ticklabels, overall x and y var labels
                if iov_yvar == 0:
                    # Column titles (overall x labels)
                    ax.set_title(ov_xvar_strnow + str(ov_xvar_val))
                    ax.axes.xaxis.set_visible(False)
                elif iov_yvar == len(ov_yvar_vals) - 1:
                    ax.set_xticks(np.arange(len(allvarvals[mat_xvar])))
                    ax.set_xticklabels(allvarvals[mat_xvar])
                else:
                    ax.axes.xaxis.set_visible(False)

                if iov_xvar == 0:
                    ax.set_yticks(np.arange(len(allvarvals[mat_yvar])))
                    ax.set_yticklabels(allvarvals[mat_yvar][::-1])
                    # Row titles (overall y labels)
                    ax.set_ylabel(ov_yvar_strnow + str(ov_yvar_val), rotation=-90, labelpad=542)
                    ax.yaxis.set_label_position("right")
                else:
                    ax.axes.yaxis.set_visible(False)

                # - Set matrix x and y var labels
                fig.text(0.5, 0.015, mat_xvar_strnow, ha='center', va='center')
                fig.text(0.085, 0.5, mat_yvar_strnow, ha='center', va='center', rotation='vertical')

                # - Add elim time annotations
                anno_textcol_threshold = ((dfed['True_Prevalence_elim_day'] / 365).min() + (dfed['True_Prevalence_elim_day'] / 365).max())/2
                for i in range(len(allvarvals[mat_yvar])):
                    for j in range(len(allvarvals[mat_xvar])):
                        if np.isnan(matnow[i, j]):
                            colornow = anno_text_colors[1]
                        else:
                            colornow = anno_text_colors[int(matnow[i, j] >= anno_textcol_threshold)]
                        annonow = round(round(matnow[i, j] / 0.1) * 0.1, 2)
                        text = ax.text(j, i, str(annonow), ha="center", va="center",
                                       color=colornow, fontsize=8)

        # - Add colorbar
        fig.subplots_adjust(right=0.9)
        cbar_ax = fig.add_axes([0.95, 0.1, 0.02, 0.8])  # left, bottom, width, height
        fig.colorbar(im, cax=cbar_ax)

        # - Save fig
        fig_file_png = os.path.join(fig_dir, wi_name + '_elim_days.png')
        fig_file_pdf = os.path.join(fig_dir, wi_name + '_elim_days.pdf')
        plt.savefig(fig_file_png, bbox_inches="tight", dpi=300)
        plt.savefig(fig_file_pdf, bbox_inches="tight", dpi=300)
        plt.show()
