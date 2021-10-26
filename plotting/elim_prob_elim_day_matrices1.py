import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rcParams
import numpy as np
import os
import pandas as pd
import plotly.colors as colors

params = {'axes.labelsize': 11,
          'axes.titlesize': 16,
          'xtick.labelsize': 11,
          'ytick.labelsize': 11}
rcParams.update(params)

rcParams['pdf.fonttype'] = 42

##
# -------- Set figures to plot
plot_elim_probs = 1
anno_ep_all_cols = 1  # 0 = don't annotate, 1 = only annotate bottom square, 2 = annotate all squares in column
anno_ep_all_rows = 0  # 0 = don't annotate, 1 = only annotate left square, 2 = annotate all squares in row

plot_elim_days = 0

outer_axes = 'rc_dord1'  # choose: 'rc_dord1', 'rr0orrr20_sneorse2'

##
# -------- Define fxns and colormaps
greens_full = cm.Greens
greens = greens_full

# - OLD PLOTLY COLORS
# greens_full = colors.get_colorscale('greens')
# greens = greens_full[1:]
# for i in range(0, len(greens)):
#     greens[i][0] = i / (len(greens) - 1)
# greens_r_full = colors.get_colorscale('greens_r')
# greens_r = greens_r_full[:-1]
# for i in range(0, len(greens_r)):
#     greens_r[i][0] = i / (len(greens_r) - 1)

anno_size_1 = 18
anno_size_2 = 16
anno_y_offset_1 = -0.26
anno_y_offset_2 = -0.3
anno_x_offset_1 = -0.35
anno_x_offset_2 = -0.3


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
    # 'spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2',
    # 'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne',
    # 'spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2',
    # 'spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne',
    # 'spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2',
    # 'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne',
    # 'spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2',
    # 'spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne',
    # 'spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2'
]
num_sweep_vars_ls = [
    4,
    # 4,
    # 4,
    # 4,
    # 4, 4, 4, 4, 4, 4
]
drive_types_ls = [
    'classic',
    # 'integral',
    # 'classic',
    # 'integral',
    # 'classic', 'integral', 'classic', 'integral', 'classic', 'integral'
]
data_dir = '..\\csvs'
fig_dir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\presentations_writeups\\' \
          'gene_drive_paper\\figures\\elim_prob_day_matrices'  # \\' + outer_axes + '_outside'
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
    if ov_xvar == 'rc':
        ov_xvar_strnow = 'Transmission blocking\n(' + ov_xvar + ') = '
    elif (ov_xvar == 'd') or (ov_xvar == 'd1'):
        ov_xvar_strnow = 'Drive efficiency\n(' + ov_xvar + ') = '
    elif (ov_xvar == 'sne') or (ov_xvar == 'se2'):
        ov_xvar_strnow = 'Fitness cost\n(' + ov_xvar + ') = '
    elif (ov_xvar == 'rr0') or (ov_xvar == 'rr20'):
        ov_xvar_strnow = 'Initial resistance\n (' + ov_xvar + ') = '
    else:
        ov_xvar_strnow = ov_xvar + ' ='

    if ov_yvar == 'rc':
        ov_yvar_strnow = 'Transmission blocking\n(' + ov_yvar + ') = '
    elif (ov_yvar == 'd') or (ov_yvar == 'd1'):
        ov_yvar_strnow = 'Drive efficiency\n(' + ov_yvar + ') = '
    elif (ov_yvar == 'sne') or (ov_yvar == 'se2'):
        ov_yvar_strnow = 'Fitness cost\n(' + ov_yvar + ') = '
    elif (ov_yvar == 'rr0') or (ov_yvar == 'rr20'):
        ov_yvar_strnow = 'Initial resistance\n (' + ov_yvar + ') = '
    else:
        ov_yvar_strnow = ov_yvar + ' ='

    if mat_xvar == 'rc':
        mat_xvar_strnow = 'Transmission blocking (' + mat_xvar + ')'
    elif (mat_xvar == 'd') or (mat_xvar == 'd1'):
        mat_xvar_strnow = 'Drive efficiency (' + mat_xvar + ')'
    elif (mat_xvar == 'sne') or (mat_xvar == 'se2'):
        mat_xvar_strnow = 'Fitness cost (' + mat_xvar + ')'
    elif (mat_xvar == 'rr0') or (mat_xvar == 'rr20'):
        mat_xvar_strnow = 'Initial resistance (' + mat_xvar + ')'
    else:
        mat_xvar_strnow = mat_xvar

    if mat_yvar == 'rc':
        mat_yvar_strnow = 'Transmission blocking (' + mat_yvar + ')'
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

        # - Initialize subplots/axes/lists
        iaxis = 1
        subplots = []
        incr_cols = []
        decr_cols = []
        wiggly_cols = []
        incr_rows = []
        decr_rows = []
        wiggly_rows = []

        dfesm = dfe[dfe[mat_xvar].isin(allvarvals[mat_xvar]) &
                    dfe[mat_yvar].isin(allvarvals[mat_yvar])]

        textcolors = ["black", "white"]
        textcolthreshold = 0.5
        valfmt = "{x:.1f}"

        fig, axes = plt.subplots(nrows=len(ov_yvar_vals), ncols=len(ov_xvar_vals),
                                 figsize=(10, 6), sharex=True, sharey=True)

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
                matnow = matnow.values
                matnow = np.flipud(matnow)

                # - Categorize columns into monotonic increasing, mt decreasing, wiggly
                incr_colsnow = []
                decr_colsnow = []
                wiggly_colsnow = []
                for icol in range(0, matnow.shape[1]):
                    colnow = matnow[:, icol]
                    incr_colsnow.append(monotonic_increasing(colnow))
                    decr_colsnow.append(monotonic_decreasing(colnow))
                    wiggly_colsnow.append(wiggly(colnow))
                incr_cols.append(incr_colsnow)
                decr_cols.append(decr_colsnow)
                wiggly_cols.append(wiggly_colsnow)

                # - Categorize rows into monotonic increasing, mt decreasing, wiggly
                incr_rowsnow = []
                decr_rowsnow = []
                wiggly_rowsnow = []
                for irow in range(0, matnow.shape[0]):
                    rownow = matnow[irow, :]
                    incr_rowsnow.append(monotonic_increasing(rownow))
                    decr_rowsnow.append(monotonic_decreasing(rownow))
                    wiggly_rowsnow.append(wiggly(rownow))
                incr_rows.append(incr_rowsnow)
                decr_rows.append(decr_rowsnow)
                wiggly_rows.append(wiggly_rowsnow)

                ax = axes[iov_yvar][iov_xvar]

                # Plot heat map
                im = ax.imshow(matnow, cmap=greens, vmin=0, vmax=1)

                # Set ticks and ticklabels
                ax.set_xticks(np.arange(len(allvarvals[mat_xvar])))
                ax.set_yticks(np.arange(len(allvarvals[mat_yvar])))
                ax.set_xticklabels(allvarvals[mat_xvar])
                ax.set_yticklabels(allvarvals[mat_yvar])
                ax.set_title(ov_xvar_strnow + str(ov_xvar_val))

                # - Row titles
                # ov_yvar_strnow + str(ov_yvar_val)

                # - Add annotations
                for i in range(len(allvarvals[mat_yvar])):
                    for j in range(len(allvarvals[mat_xvar])):
                        colornow = textcolors[int(im.norm(matnow[i, j]) > textcolthreshold)]
                        annonow = round(round(matnow[i, j]/0.05)*0.05, 2)
                        text = ax.text(j, i, str(annonow),
                                       # "{:2f}".format(annonow),
                                       ha="center", va="center", color=colornow,
                                       fontsize=8)

        # Add colorbar
        fig.subplots_adjust(right=0.9)
        cbar_ax = fig.add_axes([0.95, 0.1, 0.02, 0.8])  # left, bottom, width, height
        fig.colorbar(im, cax=cbar_ax)

        # - Update annotations for all subplots
        # symbol_anno_color = 'orange'
        # for isp, subplot in enumerate(subplots):
        #     fig.layout.annotations += subplots[isp].layout.annotations

        #     # - Annotate columns (mt increasing, mt decreasing, wiggly)
        #     mono_decr_xs = np.where(decr_cols[isp])[0]
        #     mono_incr_xs = np.where(incr_cols[isp])[0]
        #     wiggly_xs = np.where(wiggly_cols[isp])[0]
        #     if anno_ep_all_cols == 1:
        #         if len(mono_decr_xs) > 0:
        #             for anno_x in mono_decr_xs:
        #                 fig.add_annotation(x=anno_x, y=anno_y_offset_1,
        #                                    xref='x' + str(isp+1), yref='y' + str(isp+1),
        #                                    font=dict(size=anno_size_1, color=symbol_anno_color),
        #                                    text=r'$\blacktriangledown$', showarrow=False)
        #         if len(mono_incr_xs) > 0:
        #             for anno_x in mono_incr_xs:
        #                 fig.add_annotation(x=anno_x, y=anno_y_offset_1,
        #                                    xref='x' + str(isp+1), yref='y' + str(isp+1),
        #                                    font=dict(size=anno_size_1, color=symbol_anno_color),
        #                                    text=r'$\blacktriangle$', showarrow=False)
        #         if len(wiggly_xs) > 0:
        #             for anno_x in wiggly_xs:
        #                 fig.add_annotation(x=anno_x, y=anno_y_offset_1,
        #                                    xref='x' + str(isp+1), yref='y' + str(isp+1),
        #                                    font=dict(size=anno_size_1, color=symbol_anno_color),
        #                                    text='<b>~</b>', showarrow=False)
        #     elif anno_ep_all_cols == 2:
        #         if len(mono_decr_xs) > 0:
        #             for anno_x in mono_decr_xs:
        #                 for anno_y in range(0, len(allvarvals[mat_yvar])):
        #                     fig.add_annotation(x=anno_x, y=anno_y+anno_y_offset_2,
        #                                        xref='x' + str(isp+1), yref='y' + str(isp+1),
        #                                        font=dict(size=anno_size_2, color=symbol_anno_color),
        #                                        text=r'$\blacktriangledown$', showarrow=False)
        #         if len(mono_incr_xs) > 0:
        #             for anno_x in mono_incr_xs:
        #                 for anno_y in range(0, len(allvarvals[mat_yvar])):
        #                     fig.add_annotation(x=anno_x, y=anno_y+anno_y_offset_2,
        #                                        xref='x' + str(isp+1), yref='y' + str(isp+1),
        #                                        font=dict(size=anno_size_2, color=symbol_anno_color),
        #                                        text=r'$\blacktriangle$', showarrow=False)
        #         if len(wiggly_xs) > 0:
        #             for anno_x in wiggly_xs:
        #                 for anno_y in range(0, len(allvarvals[mat_yvar])):
        #                     fig.add_annotation(x=anno_x, y=anno_y+anno_y_offset_2,
        #                                        xref='x' + str(isp+1), yref='y' + str(isp+1),
        #                                        font=dict(size=anno_size_2, color=symbol_anno_color),
        #                                        text='<b>~</b>', showarrow=False)

        #     # - Annotate rows (mt increasing, mt decreasing, wiggly)
        #     mono_decr_ys = np.where(decr_rows[isp])[0]
        #     mono_incr_ys = np.where(incr_rows[isp])[0]
        #     wiggly_ys = np.where(wiggly_rows[isp])[0]
        #     if anno_ep_all_rows == 1:
        #         if len(mono_decr_ys) > 0:
        #             for anno_y in mono_decr_ys:
        #                 fig.add_annotation(x=anno_x_offset_1, y=anno_y,
        #                                    xref='x' + str(isp+1), yref='y' + str(isp+1),
        #                                    font=dict(size=anno_size_1, color=symbol_anno_color),
        #                                    text=r'$\blacktriangledown$', showarrow=False)
        #         if len(mono_incr_ys) > 0:
        #             for anno_y in mono_incr_ys:
        #                 fig.add_annotation(x=anno_x_offset_1, y=anno_y,
        #                                    xref='x' + str(isp+1), yref='y' + str(isp+1),
        #                                    font=dict(size=anno_size_1, color=symbol_anno_color),
        #                                    text=r'$\blacktriangle$', showarrow=False)
        #         if len(wiggly_ys) > 0:
        #             for anno_y in wiggly_ys:
        #                 fig.add_annotation(x=anno_x_offset_1, y=anno_y,
        #                                    xref='x' + str(isp+1), yref='y' + str(isp+1),
        #                                    font=dict(size=anno_size_1, color=symbol_anno_color),
        #                                    text='<b>~</b>', showarrow=False)
        #     elif anno_ep_all_rows == 2:
        #         if len(mono_decr_ys) > 0:
        #             for anno_y in mono_decr_ys:
        #                 for anno_x in range(0, len(allvarvals[mat_xvar])):
        #                     fig.add_annotation(y=anno_y, x=anno_x+anno_x_offset_2,
        #                                        xref='x' + str(isp+1), yref='y' + str(isp+1),
        #                                        font=dict(size=anno_size_2, color=symbol_anno_color),
        #                                        text=r'$\blacktriangledown$', showarrow=False)
        #         if len(mono_incr_ys) > 0:
        #             for anno_y in mono_incr_ys:
        #                 for anno_x in range(0, len(allvarvals[mat_xvar])):
        #                     fig.add_annotation(y=anno_y, x=anno_x+anno_x_offset_2,
        #                                        xref='x' + str(isp+1), yref='y' + str(isp+1),
        #                                        font=dict(size=anno_size_2, color=symbol_anno_color),
        #                                        text=r'$\blacktriangle$', showarrow=False)
        #         if len(wiggly_ys) > 0:
        #             for anno_y in wiggly_ys:
        #                 for anno_x in range(0, len(allvarvals[mat_xvar])):
        #                     fig.add_annotation(y=anno_y, x=anno_x+anno_x_offset_2,
        #                                        xref='x' + str(isp+1), yref='y' + str(isp+1),
        #                                        font=dict(size=anno_size_2, color=symbol_anno_color),
        #                                        text='<b>~</b>', showarrow=False)

        # fig.write_image(fig_dir + '/' + wi_name + '_elim_probs.png', width=7 * 300, height=4 * 300, scale=1)
        fig_file_png = os.path.join(fig_dir, wi_name + '_elim_probs.png')
        fig_file_pdf = os.path.join(fig_dir, wi_name + '_elim_probs.pdf')
        plt.savefig(fig_file_png, bbox_inches="tight", dpi=300)
        plt.savefig(fig_file_pdf, bbox_inches="tight", dpi=300)
        plt.show()

    ##
    # -------- Create elim day matrix
    if plot_elim_days == 1:

        # - Initialize subplots/axes
        iaxis = 1
        subplots = []

        dfedsm = dfed[dfed[mat_xvar].isin(allvarvals[mat_xvar]) &
                      dfed[mat_yvar].isin(allvarvals[mat_yvar])]

        for ov_yvar_val in ov_yvar_vals:
            for ov_xvar_val in ov_xvar_vals:

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
                matnow = (matnow / 365).round(1)  # .astype('Int64')

                # - Create annotated heatmap
                subplots.append(ff.create_annotated_heatmap(
                    z=matnow.values,
                    x=list(range(len(allvarvals[mat_xvar]))),
                    y=list(range(len(allvarvals[mat_yvar]))),
                    zmin=(dfed['True_Prevalence_elim_day'] / 365).min(),
                    zmax=(dfed['True_Prevalence_elim_day'] / 365).max(),
                    showscale=True,
                    colorscale=greens)
                )

                # - Update annotation axes
                for annot in subplots[-1]['layout']['annotations']:
                    annot['xref'] = 'x' + str(iaxis)
                    annot['yref'] = 'y' + str(iaxis)
                iaxis = iaxis + 1

        # - Set up subplot framework
        fig = make_subplots(
            rows=len(ov_yvar_vals), cols=len(ov_xvar_vals),
            shared_xaxes=True,
            shared_yaxes=True,
            column_titles=[ov_xvar_strnow + str(val) for val in ov_xvar_vals],
            row_titles=[ov_yvar_strnow + str(val) for val in ov_yvar_vals],
            x_title=mat_xvar_strnow,
            y_title=mat_yvar_strnow,
            horizontal_spacing=0.015,
            vertical_spacing=0.02
        )

        # - Change font size for row_titles/column_titles/x_title/y_title only
        for i in fig['layout']['annotations']:
            i['font'] = dict(size=24, color='black')

        # - Create each subplot
        isp = 0
        for irow, ov_yvar_val in enumerate(ov_yvar_vals):
            for icol, ov_xvar_val in enumerate(ov_xvar_vals):
                fig.add_trace(subplots[isp].data[0], row=irow + 1, col=icol + 1)
                isp = isp + 1

        # - Update annotations for all subplots
        for isp, subplot in enumerate(subplots):
            fig.layout.annotations += subplots[isp].layout.annotations

        # - Update fig layout and subplot axes
        fig.update_xaxes(
            ticklen=10,
            tickmode='array',
            tickvals=list(range(len(allvarvals[mat_xvar]))),
            ticktext=[str(val) for val in allvarvals[mat_xvar]]
        )
        fig.update_yaxes(
            ticklen=10,
            tickmode='array',
            tickvals=list(range(len(allvarvals[mat_yvar]))),
            ticktext=[str(val) for val in allvarvals[mat_yvar]]
        )
        fig.update_layout(
            # margin=dict(l=60, r=50, b=50, t=30),
            font=dict(size=20)
        )

        fig.show()
        fig.write_image(fig_dir + '/' + wi_name + '_elim_days.pdf', width=7 * 300, height=4 * 300, scale=5)
        fig.write_image(fig_dir + '/' + wi_name + '_elim_days.png', width=7 * 300, height=4 * 300, scale=5)
