##
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42

# -------- Setup params/datasets
data_dir = '..\\csvs'
base_fig_dir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\presentations_writeups\\gene_drive_paper\\figures'
num_seeds = 20  # num of seeds per sim
num_yrs = 8
released_mosqs = True
released_day = 180
itn_distrib_days = [180, 180 + 3 * 365, 180 + 6 * 365]
fc_vals = [0, 0.1, 0.2, 0.3, 0.4, 0.5]

# -- classic exps
drive_type = 'classic'
eff_allele = 'a1'
fc_var = 'sne'
# - VC and GM
distrib_itns = True
wi_name = 'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne'
# fc_type = 'increase'
# const_var_vals = {'rr0': 0.0, 'd': 1.0, 'rc': 0.5}
# const_var_vals = {'rr0': 0.0, 'd': 0.95, 'rc': 0.5}
# fc_type = 'decrease'
# const_var_vals = {'rr0': 0.0, 'd': 0.9, 'rc': 0.9}
# const_var_vals = {'rr0': 0.1, 'd': 0.9, 'rc': 0.8}
fc_type = 'incdec'
# const_var_vals = {'rr0': 0.001, 'd': 0.95, 'rc': 0.7}
# const_var_vals = {'rr0': 0.01, 'd': 0.95, 'rc': 0.7}
const_var_vals = {'rr0': 0, 'd': 0.9, 'rc': 0.6}
# - GM only
# distrib_itns = False
# -- integral exps

##
# -------- Set fig dir
fig_dir = base_fig_dir + '\\' + wi_name
os.makedirs(fig_dir, exist_ok=True)

##
# -------- Load data
dfifull = pd.read_csv(os.path.join(data_dir, 'dfi_' + wi_name + '.csv'))
dfafull = pd.read_csv(os.path.join(data_dir, 'dfa_' + wi_name + '.csv'))

##
# -------- Subset data
dfi = dfifull
dfa = dfafull
allvardefsnow = {k: v for k, v in const_var_vals.items() if k not in [fc_var]}
file_prefix = ''
for k, v in allvardefsnow.items():
    dfi = dfi[dfi[k] == v]
    dfi.drop(columns=[k], inplace=True)
    dfa = dfa[dfa[k] == v]
    dfa.drop(columns=[k], inplace=True)
    file_prefix = file_prefix + k + str(v) + '_'
file_prefix = fc_type + '_' + file_prefix
# if drive_type == 'classic':
#     dfi = dfifull[(dfifull['rc'] == const_var_vals['rc']) &
#                   (dfifull['rr0'] == const_var_vals['rr0']) &
#                   (dfifull['d'] == const_var_vals['d'])]
#     dfa = dfafull[(dfafull['rc'] == const_var_vals['rc']) &
#                   (dfafull['rr0'] == const_var_vals['rr0']) &
#                   (dfafull['d'] == const_var_vals['d'])]

##
# ---- Plot
fig, axes = plt.subplots(len(fc_vals), 1, figsize=(12, 12), sharex=True)
fig.subplots_adjust(right=0.9)

for iax, ax in enumerate(axes):
    fcnow = fc_vals[iax]
    dfinow = dfi[dfi[fc_var] == fcnow]
    dfanow = dfa[dfa[fc_var] == fcnow]

    twin1 = ax.twinx()
    twin2 = ax.twinx()
    twin3 = ax.twinx()

    # Offset the right spine of twin2.  The ticks and label have already been
    # placed on the right by twinx above.
    twin2.spines['right'].set_position(("axes", 1.07))
    twin3.spines['right'].set_position(("axes", 1.15))

    p1, = ax.plot(dfinow['Time'], dfinow['PfHRP2 Prevalence'], "k-", label="Prev")
    ax.fill_between(dfinow['Time'],
                    dfinow['PfHRP2 Prevalence'] - 1.96 * dfinow['PfHRP2 Prevalence_std'] / np.sqrt(num_seeds),
                    dfinow['PfHRP2 Prevalence'] + 1.96 * dfinow['PfHRP2 Prevalence_std'] / np.sqrt(num_seeds),
                    alpha=0.3, color='k')
    p2, = twin1.plot(dfinow['Time'], dfinow['Adult Vectors'], "r-", label="AV")
    twin1.fill_between(dfinow['Time'],
                       dfinow['Adult Vectors'] - 1.96 * dfinow['Adult Vectors_std'] / np.sqrt(num_seeds),
                       dfinow['Adult Vectors'] + 1.96 * dfinow['Adult Vectors_std'] / np.sqrt(num_seeds),
                       alpha=0.3, color='r')
    p3, = twin2.plot(dfinow['Time'], dfinow['Infectious Vectors'], "g-", label="IVF")
    twin2.fill_between(dfinow['Time'],
                       dfinow['Infectious Vectors'] - 1.96 * dfinow['Infectious Vectors_std'] / np.sqrt(num_seeds),
                       dfinow['Infectious Vectors'] + 1.96 * dfinow['Infectious Vectors_std'] / np.sqrt(num_seeds),
                       alpha=0.3, color='g')
    p4, = twin3.plot(dfanow['Time'], dfanow[eff_allele], "b-", label="EF")
    twin3.fill_between(dfanow['Time'],
                       dfanow[eff_allele] - 1.96 * dfanow[eff_allele + '_std'] / np.sqrt(num_seeds),
                       dfanow[eff_allele] + 1.96 * dfanow[eff_allele + '_std'] / np.sqrt(num_seeds),
                       alpha=0.3, color='b')

    ax.set_xlim(0, 365 * num_yrs)
    ax.set_ylim(0, 0.6)
    twin1.set_ylim(0, 4000)
    twin2.set_ylim(0, 0.1)
    twin3.set_ylim(0, 1)

    ax.set_xlabel("Time")
    ax.set_ylabel("PfHRP2 Prevalence")
    twin1.set_ylabel("Adult Vectors")
    twin2.set_ylabel("Infectious Vectors")
    twin3.set_ylabel("Effector Freq")

    ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())
    twin3.yaxis.label.set_color(p4.get_color())

    ax.set_title(fc_var + ' = ' + str(fcnow))

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
fig_file_png = os.path.join(fig_dir, file_prefix + 'prev_av_ivf_ef.png')
plt.savefig(fig_file_png, bbox_inches="tight", dpi=300)
plt.show()

