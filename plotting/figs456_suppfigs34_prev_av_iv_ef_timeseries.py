##
import pandas as pd
import numpy as np
import os
from matplotlib import rcParams
import matplotlib.pyplot as plt

params = {'axes.labelsize': 10,
          'axes.titlesize': 12,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10}
rcParams.update(params)

rcParams['pdf.fonttype'] = 42

# -------- Setup params/datasets
data_dir = '..\\csvs'
base_fig_dir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\presentations_writeups\\' \
               'gene_drive_paper\\figures\\prev_av_iv_timeseries'
num_seeds = 20  # num of seeds per sim
num_yrs = 8
released_mosqs = True
released_day = 180
itn_distrib_days = [180, 180 + 3 * 365, 180 + 6 * 365]
fig_to_plot = 'm10'  # figure number: m = main, s = supp


# ---- FITNESS COST ----

if fig_to_plot == 's3':
    drive_type = 'classic'
    eff_allele = 'a1'
    sweep_var = 'sne'
    sweep_var_ln = 'Fitness cost (sne)'
    sweep_vals = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    distrib_itns = False
    wi_name = 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne'
    sweep_type_ls = ['increase', 'decrease', 'incdec']
    const_var_vals_ls = [
        # - increase
        {'rr0': 0.001, 'd': 1.0, 'rc': 0.7},
        # - decrease
        {'rr0': 0.1, 'd': 1.0, 'rc': 1.0},
        # - incdec
        {'rr0': 0.01, 'd': 0.95, 'rc': 0.9}
    ]


# ---- PRE-EXISTING RESISTANCE ----

elif fig_to_plot == 'm6':
    drive_type = 'classic'
    eff_allele = 'a1'
    sweep_var = 'rr0'
    sweep_var_ln = 'Pre-existing resistance (rr0)'
    sweep_vals = [0.0, 0.001, 0.01, 0.1]
    distrib_itns = False
    wi_name = 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne'
    sweep_type_ls = ['decrease']
    const_var_vals_ls = [
        {'d': 0.9, 'rc': 0.9, 'sne': 0.2}
    ]


# ---- DRIVE EFFICIENCY ----

elif fig_to_plot == 'm5':
    drive_type = 'classic'
    eff_allele = 'a1'
    sweep_var = 'd'
    sweep_var_ln = 'Drive efficiency (d)'
    sweep_vals = [0.9, 0.95, 1]
    distrib_itns = False
    wi_name = 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne'
    sweep_type_ls = ['increase']
    const_var_vals_ls = [
        {'rc': 0.9, 'sne': 0.3, 'rr0': 0.01}
    ]

elif fig_to_plot == 's4':
    drive_type = 'integral'
    eff_allele = 'b1'
    sweep_var = 'd1'
    sweep_var_ln = 'Drive efficiency (d1)'
    sweep_vals = [0.9, 0.95, 1]
    distrib_itns = True
    wi_name = 'spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2'
    sweep_type_ls = ['decrease']
    const_var_vals_ls = [
        {'rc': 0.7, 'se2': 0.3, 'rr20': 0.1}
    ]


# ---- TRANSMISSION-BLOCKING EFFECTIVENESS ----

elif fig_to_plot == 'm4':
    drive_type = 'classic'
    eff_allele = 'a1'
    sweep_var = 'rc'
    sweep_var_ln = 'Transmission-blocking effectiveness (rc)'
    sweep_vals = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    distrib_itns = False
    wi_name = 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne'
    sweep_type_ls = ['increase']
    const_var_vals_ls = [
        {'d': 1, 'sne': 0.4, 'rr0': 0.01}
    ]


##
# -------- Set fig dir
fig_dir = base_fig_dir + '\\' + wi_name
os.makedirs(fig_dir, exist_ok=True)

##
# -------- Load data
dfifull = pd.read_csv(os.path.join(data_dir, 'dfi_' + wi_name + '.csv'))
dfafull = pd.read_csv(os.path.join(data_dir, 'dfa_' + wi_name + '.csv'))
dfefull = pd.read_csv(os.path.join(data_dir, 'dfe_' + wi_name + '.csv'))
dfedfull = pd.read_csv(os.path.join(data_dir, 'dfed_' + wi_name + '.csv'))

for isweep in range(0, len(sweep_type_ls)):
    const_var_vals = const_var_vals_ls[isweep]
    sweep_type = sweep_type_ls[isweep]

    ##
    # -------- Subset data
    dfi = dfifull
    dfa = dfafull
    dfe = dfefull
    dfed = dfedfull
    allvardefsnow = {k: v for k, v in const_var_vals.items() if k not in [sweep_var]}
    file_prefix = ''
    for k, v in allvardefsnow.items():
        dfi = dfi[dfi[k] == v]
        dfi.drop(columns=[k], inplace=True)
        dfa = dfa[dfa[k] == v]
        dfa.drop(columns=[k], inplace=True)
        dfe = dfe[dfe[k] == v]
        dfe.drop(columns=[k], inplace=True)
        dfed = dfed[dfed[k] == v]
        dfed.drop(columns=[k], inplace=True)
        file_prefix = file_prefix + k + str(v) + '_'
    file_prefix = sweep_type + '_' + file_prefix

    ##
    # ---- Plot
    fig, axes = plt.subplots(len(sweep_vals), 1, figsize=(7.5, 9), sharex=True)
    fig.subplots_adjust(right=0.95, hspace=0.4)

    for iax, ax in enumerate(axes):
        sweepnow = sweep_vals[iax]
        dfinow = dfi[dfi[sweep_var] == sweepnow]
        dfanow = dfa[dfa[sweep_var] == sweepnow]
        dfenow = dfe[dfe[sweep_var] == sweepnow]
        dfednow = dfed[(dfed[sweep_var] == sweepnow) & (dfed['True_Prevalence_elim'] == True)]
        epnow = dfenow['True_Prevalence_elim'].sum() / num_seeds
        ednow = dfednow['True_Prevalence_elim_day'].mean()

        twin1 = ax.twinx()
        twin2 = ax.twinx()
        twin3 = ax.twinx()
        twin4 = ax.twinx()

        # Offset the right spine of twin2.  The ticks and label have already been
        # placed on the right by twinx above.
        twin2.spines['right'].set_position(("axes", 1.105))
        twin3.spines['right'].set_position(("axes", 1.22))
        twin4.spines['right'].set_position(("axes", 1.315))

        p1, = ax.plot(dfinow['Time'], dfinow['True Prevalence'], color='k', label='Prev')
        ax.fill_between(dfinow['Time'],
                        dfinow['True Prevalence'] - 1.96 * dfinow['True Prevalence_std'] / np.sqrt(num_seeds),
                        dfinow['True Prevalence'] + 1.96 * dfinow['True Prevalence_std'] / np.sqrt(num_seeds),
                        alpha=0.3, color='k')
        p2, = twin1.plot(dfinow['Time'], dfinow['Adult Vectors'], color='tab:orange', label="AV")
        twin1.fill_between(dfinow['Time'],
                           dfinow['Adult Vectors'] - 1.96 * dfinow['Adult Vectors_std'] / np.sqrt(num_seeds),
                           dfinow['Adult Vectors'] + 1.96 * dfinow['Adult Vectors_std'] / np.sqrt(num_seeds),
                           alpha=0.3, color='tab:orange')
        p3, = twin2.plot(dfinow['Time'], dfinow['Infectious Vectors'], color='tab:green', label="IVF")
        twin2.fill_between(dfinow['Time'],
                           dfinow['Infectious Vectors'] - 1.96 * dfinow['Infectious Vectors_std'] / np.sqrt(num_seeds),
                           dfinow['Infectious Vectors'] + 1.96 * dfinow['Infectious Vectors_std'] / np.sqrt(num_seeds),
                           alpha=0.3, color='tab:green')
        p4, = twin3.plot(dfinow['Time'], dfinow['Infectious Vectors']*dfinow['Adult Vectors'], color='tab:red', label="IVN")
        p5, = twin4.plot(dfanow['Time'], dfanow[eff_allele], color='tab:blue', label="EF")
        twin4.fill_between(dfanow['Time'],
                           dfanow[eff_allele] - 1.96 * dfanow[eff_allele + '_std'] / np.sqrt(num_seeds),
                           dfanow[eff_allele] + 1.96 * dfanow[eff_allele + '_std'] / np.sqrt(num_seeds),
                           alpha=0.3, color='tab:blue')

        ax.set_xlim(0, 365 * num_yrs)
        ax.set_ylim(0, 1)
        twin1.set_ylim(0, 4000)
        twin2.set_ylim(0, 0.1)
        twin3.set_ylim(0, 100)
        twin4.set_ylim(0, 1)

        ax.set_ylabel("Prevalence")
        twin1.set_ylabel("Adult Vectors")
        twin2.set_ylabel("Inf. Vector Frac.")
        twin3.set_ylabel("Inf. Vectors")
        twin4.set_ylabel("Effector Freq")

        ax.yaxis.label.set_color(p1.get_color())
        twin1.yaxis.label.set_color(p2.get_color())
        twin2.yaxis.label.set_color(p3.get_color())
        twin3.yaxis.label.set_color(p4.get_color())
        twin4.yaxis.label.set_color(p5.get_color())

        ax.set_title(sweep_var_ln + ' = ' + str(sweepnow) +
                     ', ' + 'e.p. = ' + str(epnow) +
                     ', ' + 'e.d. = ' + "{:.1f}".format(ednow))

        tkw = dict(size=2.5, width=1.5)
        ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
        twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
        twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
        twin3.tick_params(axis='y', colors=p4.get_color(), **tkw)
        twin4.tick_params(axis='y', colors=p5.get_color(), **tkw)
        ax.tick_params(axis='x', **tkw)

        # ax.legend(handles=[p1, p2, p3, p4])

        if released_mosqs == True:
            ax.axvline(x=released_day, color='k', linestyle='--')

        if distrib_itns == True:
            for itn_day in itn_distrib_days:
                ax.axvline(x=itn_day, color='gray', linestyle=':')

    fig.text(0.5, 0.07, "Time (days)", ha='center', va='center')

    fig_file_png = os.path.join(fig_dir, file_prefix + 'prev_av_ivf_ef.png')
    fig_file_pdf = os.path.join(fig_dir, file_prefix + 'prev_av_ivf_ef.pdf')
    plt.savefig(fig_file_png, bbox_inches="tight", dpi=300)
    plt.savefig(fig_file_pdf, bbox_inches="tight", dpi=300)
    plt.show()
