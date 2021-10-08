import numpy as np
import os
import pandas as pd
import plotly.colors as colors
import plotly.figure_factory as ff
from plotly.subplots import make_subplots

greens_full = colors.get_colorscale('greens')
greens = greens_full[1:]
for i in range(0, len(greens)):
    greens[i][0] = i / (len(greens) - 1)


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


plot_elim_probs = 1
plot_elim_days = 0

# -------- Setup params/datasets
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
    4,  # 4, 4, 4,
    # 4, 4, 4, 4, 4, 4
]
drive_types_ls = [
    'classic',  # 'integral', 'classic', 'integral',
    # 'classic', 'integral', 'classic', 'integral', 'classic', 'integral'
]
data_dir = '..\\csvs'
fig_dir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\presentations_writeups\\gene_drive_paper\\figures\\elim_prob_day_matrices'
num_seeds = 20  # num of seeds per sim

for iwi, wi_name in enumerate(wi_names_ls):
    drive_type = drive_types_ls[iwi]
    num_sweep_vars = num_sweep_vars_ls[iwi]

    if num_sweep_vars == 6:
        if drive_type == 'classic':
            allvardefs = {'rc (phenotypic effectiveness)': 1, 'd (drive efficiency)': 1,
                          'rr0 (initial resistance)': 0, 'sne (fitness cost)': 0,
                          'rd': 180, 'nn': 6}
            allvarvals = {
                'rc (phenotypic effectiveness)': [0.5, 0.6, 0.7, 0.8, 0.9, 1],
                'd (drive efficiency)': [0.9, 0.95, 1],
                'rr0 (initial resistance)': [0, 0.001, 0.01, 0.1],
                'sne (fitness cost)': [0, 0.1, 0.2, 0.3, 0.4, 0.5],
                'rd': [180, 240, 300, 360, 420, 480, 545],
                'nn': [6, 12]
            }
    elif num_sweep_vars == 4:
        if drive_type == 'classic':
            allvardefs = {'rc (phenotypic effectiveness)': 1, 'd (drive efficiency)': 1,
                          'sne (fitness cost)': 0, 'rr0 (initial resistance)': 0}
            allvarvals = {
                'rc (phenotypic effectiveness)': [0.5, 0.6, 0.7, 0.8, 0.9, 1],
                'd (drive efficiency)': [0.9, 0.95, 1],
                'rr0 (initial resistance)': [0, 0.001, 0.01, 0.1],
                'sne (fitness cost)': [0, 0.1, 0.2, 0.3, 0.4, 0.5]
            }
        elif drive_type == 'integral':
            allvardefs = {'rc (phenotypic effectiveness)': 1, 'd1 (drive efficiency)': 1,
                          'se2 (fitness cost)': 0, 'rr20 (initial resistance)': 0}
            allvarvals = {
                'rc (phenotypic effectiveness)': [0.5, 0.6, 0.7, 0.8, 0.9, 1],
                'd1 (drive efficiency)': [0.9, 0.95, 1],
                'rr20 (initial resistance)': [0, 0.001, 0.01, 0.1],
                'se2 (fitness cost)': [0, 0.1, 0.2, 0.3, 0.4, 0.5]
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

    if drive_type == 'classic':
        dfe.rename(columns={'rc': 'rc (phenotypic effectiveness)',
                            'd': 'd (drive efficiency)',
                            'rr0': 'rr0 (initial resistance)',
                            'sne': 'sne (fitness cost)',
                            }, inplace=True)
        dfed.rename(columns={'rc': 'rc (phenotypic effectiveness)',
                             'd': 'd (drive efficiency)',
                             'rr0': 'rr0 (initial resistance)',
                             'sne': 'sne (fitness cost)',
                             }, inplace=True)
    elif drive_type == 'integral':
        dfe.rename(columns={'rc': 'rc (phenotypic effectiveness)',
                            'd1': 'd1 (drive efficiency)',
                            'rr20': 'rr20 (initial resistance)',
                            'se2': 'se2 (fitness cost)',
                            }, inplace=True)
        dfed.rename(columns={'rc': 'rc (phenotypic effectiveness)',
                             'd1': 'd1 (drive efficiency)',
                             'rr20': 'rr20 (initial resistance)',
                             'se2': 'se2 (fitness cost)',
                             }, inplace=True)

    ##
    # -------- Set up x/y axes
    # - Set matrix x/y vars, overall x/y vars
    if drive_type == 'classic':
        # - rc/d on outer axes
        mat_xvar = 'rr0 (initial resistance)'
        mat_yvar = 'sne (fitness cost)'
        ov_xvar = 'rc (phenotypic effectiveness)'
        ov_yvar = 'd (drive efficiency)'
        # - rr0/sne on outer axes
        # mat_xvar = 'rc (phenotypic effectiveness)'
        # mat_yvar = 'd (drive efficiency)'
        # ov_xvar = 'sne (fitness cost)'
        # ov_yvar = 'rr0 (initial resistance)'
    elif drive_type == 'integral':
        # - rc/d1 on outer axes
        mat_xvar = 'rr20 (initial resistance)'
        mat_yvar = 'se2 (fitness cost)'
        ov_xvar = 'rc (phenotypic effectiveness)'
        ov_yvar = 'd1 (drive efficiency)'
        # - rr0/sne on outer axes
        # mat_xvar = 'rc (phenotypic effectiveness)'
        # mat_yvar = 'd1 (drive efficiency)'
        # ov_xvar = 'se2 (fitness cost)'
        # ov_yvar = 'rr20 (initial resistance)'

    ov_xvar_vals = allvarvals[ov_xvar]
    ov_yvar_vals = allvarvals[ov_yvar]

    # -------- Create elim prob matrix
    if plot_elim_probs == 1:

        # - Initialize subplots/axes/lists
        iaxis = 1
        subplots = []
        incr_cols = []
        decr_cols = []
        wiggly_cols = []

        dfesm = dfe[dfe[mat_xvar].isin(allvarvals[mat_xvar]) &
                    dfe[mat_yvar].isin(allvarvals[mat_yvar])]

        for ov_yvar_val in ov_yvar_vals:
            for ov_xvar_val in ov_xvar_vals:

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

                # - try something - Compute monotonic increasing columns
                matnow = matnow.values
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

                # - Create annotated heatmap
                subplots.append(ff.create_annotated_heatmap(
                    z=matnow,
                    x=list(range(len(allvarvals[mat_xvar]))),
                    y=list(range(len(allvarvals[mat_yvar]))),
                    zmin=0,
                    zmax=1,
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
            column_titles=[ov_xvar + '=' + str(val) for val in ov_xvar_vals],
            row_titles=[ov_yvar + '=' + str(val) for val in ov_yvar_vals],
            x_title=mat_xvar,
            y_title=mat_yvar,
            horizontal_spacing=0.03,
            vertical_spacing=0.03
        )

        # - Create each subplot
        isp = 0
        for irow, ov_yvar_val in enumerate(ov_yvar_vals):
            for icol, ov_xvar_val in enumerate(ov_xvar_vals):
                fig.add_trace(subplots[isp].data[0], row=irow + 1, col=icol + 1)
                isp = isp + 1

        # - Update annotations for all subplots
        for isp, subplot in enumerate(subplots):
            fig.layout.annotations += subplots[isp].layout.annotations
            # - try something - Annotate monotonically increasing/decreasing columns
            mono_decr_xs = np.where(decr_cols[isp])[0]
            mono_incr_xs = np.where(incr_cols[isp])[0]
            wiggly_xs = np.where(wiggly_cols[isp])[0]
            if len(mono_decr_xs) > 0:
                for anno_x in mono_decr_xs:
                    for anno_y in range(0, len(allvarvals[mat_yvar])):
                        fig.add_annotation(x=anno_x, y=anno_y-0.3,
                        # fig.add_annotation(x=anno_x, y=-0.3,
                                           xref='x' + str(isp+1), yref='y' + str(isp+1),
                                           font=dict(size=14, color='red'),
                                           # font=dict(size=16, color='red'),
                                           text=r'$\triangledown$', showarrow=False)
            if len(mono_incr_xs) > 0:
                for anno_x in mono_incr_xs:
                    for anno_y in range(0, len(allvarvals[mat_yvar])):
                        fig.add_annotation(x=anno_x, y=anno_y-0.3,
                        # fig.add_annotation(x=anno_x, y=-0.3,
                                           xref='x' + str(isp+1), yref='y' + str(isp+1),
                                           font=dict(size=11, color='red'),
                                           # font=dict(size=12, color='red'),
                                           text=r'$\triangle$', showarrow=False)
            if len(wiggly_xs) > 0:
                for anno_x in wiggly_xs:
                    for anno_y in range(0, len(allvarvals[mat_yvar])):
                        fig.add_annotation(x=anno_x, y=anno_y-0.3,
                        # fig.add_annotation(x=anno_x, y=-0.3,
                                           xref='x' + str(isp+1), yref='y' + str(isp+1),
                                           font=dict(size=14, color='red'),
                                           # font=dict(size=16, color='red'),
                                           text='~', showarrow=False)

        # - Update fig layout and subplot axes
        fig.update_xaxes(
            tickmode='array',
            tickvals=list(range(len(allvarvals[mat_xvar]))),
            ticktext=[str(val) for val in allvarvals[mat_xvar]]
        )
        fig.update_yaxes(
            tickmode='array',
            tickvals=list(range(len(allvarvals[mat_yvar]))),
            ticktext=[str(val) for val in allvarvals[mat_yvar]]
        )
        fig.update_layout(margin=dict(l=60, r=50, b=50, t=30))

        fig.show()
        fig.write_image(fig_dir + '/' + wi_name + '_elim_probs.pdf', width=7 * 300, height=4 * 300, scale=5)
        fig.write_image(fig_dir + '/' + wi_name + '_elim_probs.png', width=7 * 300, height=4 * 300, scale=5)

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
            column_titles=[ov_xvar + '=' + str(val) for val in ov_xvar_vals],
            row_titles=[ov_yvar + '=' + str(val) for val in ov_yvar_vals],
            x_title=mat_xvar,
            y_title=mat_yvar,
            horizontal_spacing=0.03,
            vertical_spacing=0.03
        )

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
        fig.update_layout(margin=dict(l=60, r=50, b=50, t=30))

        fig.show()
        fig.write_image(fig_dir + '/' + wi_name + '_elim_days.pdf', width=7 * 300, height=4 * 300, scale=5)
        fig.write_image(fig_dir + '/' + wi_name + '_elim_days.png', width=7 * 300, height=4 * 300, scale=5)
