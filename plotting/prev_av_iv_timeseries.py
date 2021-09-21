##
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42

# -------- Setup params/datasets
data_dir = '..\\csvs'
fig_dir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\presentations_writeups\\gene_drive_paper\\figures'
num_seeds = 20  # num of seeds per sim
num_yrs = 8
released_mosqs = True
released_day = 180
itn_distrib_days = [180, 180 + 3 * 365, 180 + 6 * 365]
fc_vals = [0, 0.1, 0.2, 0.3, 0.4, 0.5]

# -- classic exps
eff_allele = 'a1'
fc_var = 'sne'
# - VC and GM
distrib_itns = True
wi_name = 'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne'
const_var_vals = {'rr0': 0, 'd': 1, 'rc': 0.5, 'rd': 180, 'nn': 6}
# - GM only
# distrib_itns = False

# - integral exps

##
# -------- Load data
dfi = pd.read_csv(os.path.join(data_dir, 'dfi_' + wi_name + '.csv'))
dfa = pd.read_csv(os.path.join(data_dir, 'dfa_' + wi_name + '.csv'))

##
# -------- Subset data
dfi = dfi[]

allvardefsnow = {k: v for k, v in const_var_vals.items() if k not in [fc_var]}
for k, v in allvardefsnow.items():
    dfenow = dfesm[dfesm[k] == v]
    dfenow.drop(columns=[k], inplace=True)

##
# ---- Plot all in one
plotnow = 1

if plotnow == 1:
    fig, axes = plt.subplots(len(fc_vals), 1, figsize=(14, 10))
    fig.subplots_adjust(right=0.9)

    for iax, ax in enumerate(axes):
        fcnow = fc_vals[iax]

        twin1 = ax.twinx()
        twin2 = ax.twinx()
        twin3 = ax.twinx()

        # Offset the right spine of twin2.  The ticks and label have already been
        # placed on the right by twinx above.
        twin2.spines['right'].set_position(("axes", 1.07))
        twin3.spines['right'].set_position(("axes", 1.15))

        p1, = ax.plot(dfi['Time'], dfi['PfHRP2 Prevalence'], "k-", label="Prev")
        ax.fill_between(dfi['Time'],
                        dfi['PfHRP2 Prevalence'] - 1.96 * dfi['PfHRP2 Prevalence_std'] / np.sqrt(num_seeds),
                        dfi['PfHRP2 Prevalence'] + 1.96 * dfi['PfHRP2 Prevalence_std'] / np.sqrt(num_seeds),
                        alpha=0.3, color='k')
        p2, = twin1.plot(dfi['Time'], dfi['Adult Vectors'], "r-", label="AV")
        twin1.fill_between(dfi['Time'],
                           dfi['Adult Vectors'] - 1.96 * dfi['Adult Vectors_std'] / np.sqrt(num_seeds),
                           dfi['Adult Vectors'] + 1.96 * dfi['Adult Vectors_std'] / np.sqrt(num_seeds),
                           alpha=0.3, color='r')
        p3, = twin2.plot(dfi['Time'], dfi['Infectious Vectors'], "g-", label="IVF")
        twin2.fill_between(dfi['Time'],
                           dfi['Infectious Vectors'] - 1.96 * dfi['Infectious Vectors_std'] / np.sqrt(num_seeds),
                           dfi['Infectious Vectors'] + 1.96 * dfi['Infectious Vectors_std'] / np.sqrt(num_seeds),
                           alpha=0.3, color='g')
        p4, = twin3.plot(dfa['Time'], dfa[eff_allele], "b-", label="EF")
        twin3.fill_between(dfa['Time'],
                           dfa[eff_allele] - 1.96 * dfa[eff_allele + '_std'] / np.sqrt(num_seeds),
                           dfa[eff_allele] + 1.96 * dfa[eff_allele + '_std'] / np.sqrt(num_seeds),
                           alpha=0.3, color='b')

        ax.set_xlim(0, 365 * num_yrs)
        ax.set_ylim(0, 0.6)
        twin1.set_ylim(0, 4000)
        twin2.set_ylim(0, 0.1)
        twin3.set_ylim(0, 1)

        ax.set_xlabel("Time")
        ax.set_ylabel("PfHRP2 Prevalence")
        twin1.set_ylabel("Adult Vectors")
        twin2.set_ylabel("Infectious Vectors Fraction")
        twin3.set_ylabel("Effector Freq")

        ax.yaxis.label.set_color(p1.get_color())
        twin1.yaxis.label.set_color(p2.get_color())
        twin2.yaxis.label.set_color(p3.get_color())
        twin3.yaxis.label.set_color(p4.get_color())

        tkw = dict(size=4, width=1.5)
        ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
        twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
        twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
        twin3.tick_params(axis='y', colors=p4.get_color(), **tkw)
        ax.tick_params(axis='x', **tkw)

        # ax.legend(handles=[p1, p2, p3, p4])

        if released_mosqs == True:
            ax.axvline(x=released_day, color='k', linestyle='--')  # mosq release date

        if distrib_itns == True:
            for itn_day in itn_distrib_days:
                ax.axvline(x=itn_day, color='gray', linestyle=':')

        fig.tight_layout()
        fig_file_png = os.path.join(fig_dir, 'prev_av_ivf_ef_seed.png')
        plt.savefig(fig_file_png, bbox_inches="tight", dpi=300)
        plt.show()
