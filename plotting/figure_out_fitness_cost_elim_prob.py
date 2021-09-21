##
from createSimDirectoryMapBR import createSimDirectoryMap
import json
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42

##
# ------ Specify exp and params
num_yrs = 8
num_seeds = 20
released_mosqs = True
released_day = 180
distrib_itns = True
itn_distrib_days = [180, 180 + 3 * 365, 180 + 6 * 365]
drive_type = "classic"
alleles = ['a0', 'a1', 'a2']

# -- spatial, classic, GM only, EIR = 10
# wi_name = 'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_release_node_num'
# simmap = createSimDirectoryMap('0d04ced1-78f1-eb11-a9ed-b88303911bc1')
# - elim prob up w/ fitness cost
# wi_name_fig_dir = 'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_release_node_num\\epinc\\sne0pt2'
# wi_name_fig_dir = 'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_release_node_num\\epinc\\sne0pt05'
# simmapnow = simmap[(simmap['rr0'] == 0) & (simmap['d'] == 1) & (simmap['rc'] == 0.7)
#                    & (simmap['release_day'] == 180) & (simmap['num_nodes'] == 6)
#                    & (simmap['sne'] == 0.2)]
#                    # & (simmap['sne'] == 0.05)]
# - elim prob down w/ fitness cost
# wi_name_fig_dir = 'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_release_node_num\\epdec\\sne0pt2'
# wi_name_fig_dir = 'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_release_node_num\\epdec\\sne0'
# simmapnow = simmap[(simmap['rr0'] == 0.1) & (simmap['d'] == 0.9) & (simmap['rc'] == 0.9)
#                    & (simmap['release_day'] == 180) & (simmap['num_nodes'] == 6)
#                    # & (simmap['sne'] == 0.2)]
#                    & (simmap['sne'] == 0)]
# - elim prob up then down w/ fitness cost
# wi_name_fig_dir = 'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_release_node_num\\epincdec\\sne0pt2'
# wi_name_fig_dir = 'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_release_node_num\\epincdec\\sne0pt1'
# wi_name_fig_dir = 'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_release_node_num\\epincdec\\sne0'
# simmapnow = simmap[(simmap['rr0'] == 0) & (simmap['d'] == 0.95) & (simmap['rc'] == 0.8)
#                    & (simmap['release_day'] == 180) & (simmap['num_nodes'] == 6)
#                    & (simmap['sne'] == 0.2)]
#                    # & (simmap['sne'] == 0.1)]
#                    # & (simmap['sne'] == 0)]

# -- spatial, classic, GM only, EIR = 30
# - elim prob up w/ fitness cost
# wi_name = 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newrr0'
# wi_name_fig_dir = 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newrr0\\sne0pt2'
# # wi_name_fig_dir = 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newrr0\\sne0'
# simmap = createSimDirectoryMap('5e309994-860a-ec11-a9ed-b88303911bc1')
# simmapnow = simmap[(simmap['rr0'] == 0.001) & (simmap['d'] == 1) & (simmap['rc'] == 0.8)
#                    & (simmap['sne'] == 0.2)]
#                    # & (simmap['sne'] == 0)]
# - elim prob down w/ fitness cost
wi_name = 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_sne_newsne'
wi_name_fig_dir = 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_sne_newsne\\sne0pt25'
# wi_name_fig_dir = 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_sne_newsne\\sne0pt5'
simmap = createSimDirectoryMap('c7809967-c40a-ec11-a9ed-b88303911bc1')
simmapnow = simmap[(simmap['rr0'] == 0.01) & (simmap['d'] == 0.95) & (simmap['rc'] == 0.9)
                   # & (simmap['sne'] == 0.5)]
                   & (simmap['sne'] == 0.25)]

# -- spatial, classic, VC and GM, EIR = 30
# - elim prob up w/ fitness cost
# wi_name = 'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newrr0'
# # wi_name_fig_dir = 'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newrr0\\sne0pt2'
# wi_name_fig_dir = 'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newrr0\\sne0pt1'
# # wi_name_fig_dir = 'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newrr0\\sne0'
# simmap = createSimDirectoryMap('c9d8e922-860a-ec11-a9ed-b88303911bc1')
# simmapnow = simmap[(simmap['rr0'] == 0.01) & (simmap['d'] == 1) & (simmap['rc'] == 0.6)
#                    # & (simmap['sne'] == 0.2)]
#                    & (simmap['sne'] == 0.1)]
#                    # & (simmap['sne'] == 0)]

##
# ------ Load data
fpaths = simmapnow['outpath'].values

adult_vectors = pd.DataFrame()
for iseed in range(0, num_seeds):
    ic_json = json.loads(open(os.path.join(fpaths[iseed], "output", "InsetChart.json")).read())
    adult_vectors[iseed] = ic_json["Channels"]["Adult Vectors"]["Data"]
adult_vectors.reset_index(inplace=True)
adult_vectors.rename(columns={'index': 'time'}, inplace=True)

inf_vectors = pd.DataFrame()
for iseed in range(0, num_seeds):
    ic_json = json.loads(open(os.path.join(fpaths[iseed], "output", "InsetChart.json")).read())
    inf_vectors[iseed] = ic_json["Channels"]["Infectious Vectors"]["Data"]
inf_vectors.reset_index(inplace=True)
inf_vectors.rename(columns={'index': 'time'}, inplace=True)

rdt_prev = pd.DataFrame()
for iseed in range(0, num_seeds):
    ic_json = json.loads(open(os.path.join(fpaths[iseed], "output", "InsetChart.json")).read())
    rdt_prev[iseed] = ic_json["Channels"]["PfHRP2 Prevalence"]["Data"]
rdt_prev.reset_index(inplace=True)
rdt_prev.rename(columns={'index': 'time'}, inplace=True)

true_prev = pd.DataFrame()
for iseed in range(0, num_seeds):
    ic_json = json.loads(open(os.path.join(fpaths[iseed], "output", "InsetChart.json")).read())
    true_prev[iseed] = ic_json["Channels"]["True Prevalence"]["Data"]
true_prev.reset_index(inplace=True)
true_prev.rename(columns={'index': 'time'}, inplace=True)

allele_freq = pd.DataFrame()
for iseed in range(0, num_seeds):
    datatemp = pd.read_csv(os.path.join(fpaths[iseed], "output", "ReportVectorGenetics_gambiae_Female_ALLELE_FREQ.csv"))
    datatemp = datatemp[datatemp['Alleles'] != 'Y']
    datatemp = datatemp.pivot_table('VectorPopulation', ['Time', 'NodeID'], 'Alleles').reset_index()
    datatemp = datatemp.groupby('Time').sum().reset_index()
    for allele in alleles:
        allele_freq[allele + '_' + str(iseed)] = datatemp[allele] / datatemp['X']
allele_freq.reset_index(inplace=True)
allele_freq.rename(columns={'index': 'time'}, inplace=True)

##
# ------ Set fig paths
fig_dir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\model_figures\\COMPS_figures\\' + wi_name_fig_dir
os.makedirs(fig_dir, exist_ok=True)

##
# ---- Specify colors and linestyles

linestyles = ['-', '--', ':']

color_start, color_stop = 0.7, 0.1
cmap = mpl.cm.get_cmap('afmhot')

af_color_start, af_color_stop = 0.2, 0.9
af_cmap = mpl.cm.get_cmap('plasma')


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

##
# ---- Plot all in one
plotnow = 1

if plotnow == 1:
    for iseed in range(0, num_seeds):
        fig, ax = plt.subplots(1, 1, figsize=(14, 4))
        fig.subplots_adjust(right=0.9)

        twin1 = ax.twinx()
        twin2 = ax.twinx()
        twin3 = ax.twinx()

        # Offset the right spine of twin2.  The ticks and label have already been
        # placed on the right by twinx above.
        twin2.spines['right'].set_position(("axes", 1.07))
        twin3.spines['right'].set_position(("axes", 1.15))

        p1, = ax.plot(rdt_prev['time'], rdt_prev[iseed], "k-", label="Prev")
        p2, = twin1.plot(adult_vectors['time'], adult_vectors[iseed], "r-", label="AV")
        p3, = twin2.plot(inf_vectors['time'], inf_vectors[iseed], "g-", label="IVF")
        p4, = twin3.plot(allele_freq['time'], allele_freq['a1_' + str(iseed)], "b-", label="EF")

        ax.set_xlim(0, 365*num_yrs)
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
        fig_file_png = os.path.join(fig_dir, 'prev_av_ivf_ef_seed' + str(iseed) + '.png')
        plt.savefig(fig_file_png, bbox_inches="tight", dpi=300)
        plt.show()
