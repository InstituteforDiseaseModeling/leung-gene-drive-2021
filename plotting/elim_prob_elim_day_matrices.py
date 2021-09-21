import numpy as np
import os
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff
from plotly.subplots import make_subplots

plot_elim_probs = 0
plot_elim_days = 1

# -------- Setup params/datasets
wi_names_ls = [
    'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne',
    'spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2',
    'spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2',
    'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne',
    'spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne',
    'spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2',
    'spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2',
    'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne',
    'spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne',
    'spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2'
]
num_sweep_vars_ls = [
    4, 4, 4, 6,
    4, 4, 4, 6, 6, 4
    ]
drive_types_ls = [
    'classic', 'integral', 'integral', 'classic',
    'classic', 'integral', 'integral', 'classic', 'classic', 'integral'
]
data_dir = '..\\csvs'
fig_dir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\presentations_writeups\\gene_drive_paper\\figures'
num_seeds = 20  # num of seeds per sim

for iwi, wi_name in enumerate(wi_names_ls):
    drive_type = drive_types_ls[iwi]
    num_sweep_vars = num_sweep_vars_ls[iwi]

    if num_sweep_vars == 6:
        if drive_type == 'classic':
            allvardefs = {'rc': 1, 'd': 1, 'rr0': 0, 'sne': 0,
                          'rd': 180, 'nn': 6}
            allvarvals = {'rc': [1, 0.9, 0.8, 0.7, 0.6, 0.5],
                          'd': [1, 0.95, 0.9],
                          'rr0': [0, 0.001, 0.01, 0.1],
                          'sne': [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5],
                          # 'sne': [0, 0.1, 0.2, 0.3, 0.4, 0.5],
                          'rd': [180, 240, 300, 360, 420, 480, 545],
                          'nn': [6, 12]}
    elif num_sweep_vars == 4:
        if drive_type == 'classic':
            allvardefs = {'rc': 1, 'd': 1, 'sne': 0, 'rr0': 0}
            allvarvals = {'rc': [1, 0.9, 0.8, 0.7, 0.6, 0.5],
                          'd': [1, 0.95, 0.9],
                          'rr0': [0, 0.001, 0.01, 0.1],
                          'sne': [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]}
                          # 'sne': [0, 0.1, 0.2, 0.3, 0.4, 0.5]}
        elif drive_type == 'integral':
            allvardefs = {'rc': 1, 'd1': 1, 'se2': 0, 'rr20': 0}
            allvarvals = {'rc': [1, 0.9, 0.8, 0.7, 0.6, 0.5],
                          'd1': [1, 0.95, 0.9],
                          'rr20': [0, 0.001, 0.01, 0.1],
                          'se2': [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]}
                          # 'se2': [0, 0.1, 0.2, 0.3, 0.4, 0.5]}

    ##
    # -------- Load data
    # dfi = pd.read_csv(os.path.join(data_dir, 'dfi_' + wi_name + '.csv'))
    # dfa = pd.read_csv(os.path.join(data_dir, 'dfa_' + wi_name + '.csv'))
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
        mat_xvar = 'rr0'
        mat_yvar = 'sne'
        ov_xvar = 'rc'
        ov_yvar = 'd'
    elif drive_type == 'integral':
        mat_xvar = 'rr20'
        mat_yvar = 'se2'
        ov_xvar = 'rc'
        ov_yvar = 'd1'

    ov_xvar_vals = allvarvals[ov_xvar]
    ov_yvar_vals = allvarvals[ov_yvar]

    # -------- Create elim prob matrix
    if plot_elim_probs == 1:

        # - Initialize subplots/axes
        iaxis = 1
        subplots = []

        dfesm = dfe[dfe[mat_xvar].isin(allvarvals[mat_xvar]) &
                    dfe[mat_yvar].isin(allvarvals[mat_yvar])]

        for ov_yvar_val in ov_yvar_vals:
            for ov_xvar_val in ov_xvar_vals:

                # - Compute heatmap values
                dfenow = dfesm
                allvardefsnow = {k: v for k, v in allvardefs.items() if k not in [mat_xvar, mat_yvar, ov_xvar, ov_yvar]}
                if len(allvardefsnow) > 0:
                    for k, v in allvardefsnow.items():
                        dfenow = dfenow[dfenow[k] == v]
                        dfenow.drop(columns=[k], inplace=True)
                dfenow = dfenow[dfenow[ov_xvar] == ov_xvar_val]
                dfenow = dfenow[dfenow[ov_yvar] == ov_yvar_val]
                dfenow.drop(columns=[ov_xvar, ov_yvar], inplace=True)
                dfenownow = (dfenow.groupby([mat_xvar, mat_yvar])['True_Prevalence_elim'].sum() / num_seeds).reset_index()
                matnow = dfenownow.pivot_table(index=[mat_yvar], columns=[mat_xvar], values='True_Prevalence_elim')

                # - Create annotated heatmap
                subplots.append(ff.create_annotated_heatmap(
                    z=matnow.values,
                    x=list(range(len(allvarvals[mat_xvar]))),
                    y=list(range(len(allvarvals[mat_yvar]))),
                    zmin=0,
                    zmax=1,
                    showscale=True,
                    colorscale='YlOrBr_r')
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
        fig.write_image(fig_dir + '/' + wi_name + '_elim_probs.pdf', width=7*300, height=4*300, scale=5)
        fig.write_image(fig_dir + '/' + wi_name + '_elim_probs.png', width=7*300, height=4*300, scale=5)

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
                dfednow = dfedsm
                allvardefsnow = {k: v for k, v in allvardefs.items() if k not in [mat_xvar, mat_yvar, ov_xvar, ov_yvar]}
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
                    colorscale='YlOrBr')
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
        fig.write_image(fig_dir + '/' + wi_name + '_elim_days.pdf', width=7*300, height=4*300, scale=5)
        fig.write_image(fig_dir + '/' + wi_name + '_elim_days.png', width=7*300, height=4*300, scale=5)
