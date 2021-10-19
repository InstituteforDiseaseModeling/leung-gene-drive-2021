##
from createSimDirectoryMapBR import createSimDirectoryMap
import json
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import pandas as pd

params = {'legend.fontsize': 16,
          'axes.labelsize': 20,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18}
pylab.rcParams.update(params)

rcParams['pdf.fonttype'] = 42

##
# ------ Specify exp and params
released_day = 180
eir_labels = ['EIR = 10 (low transmission)', 'EIR = 30 (moderate transmission)', 'EIR = 80 (high transmission)']
simmaps = [[]] * 3
num_seeds = 16
# - spatial, integral, EIR = 10 serialization
simmaps[0] = createSimDirectoryMap('2b699ff7-f80a-ec11-a9ed-b88303911bc1')
# - spatial, integral, EIR = 30 serialization
simmaps[1] = createSimDirectoryMap('f9fbc98e-6821-ec11-9ecd-9440c9bee941')
# - spatial, integral, EIR = 80 serialization
simmaps[2] = createSimDirectoryMap('bf2fcc45-f90a-ec11-a9ed-b88303911bc1')

##
# ------ Load data
adult_vector_means = [[]] * 3
for isim, simmap in enumerate(simmaps):
    fpaths = simmap['outpath'].values

    adult_vectors = pd.DataFrame()
    for iseed in range(0, num_seeds):
        ic_json = json.loads(open(os.path.join(fpaths[iseed], "output", "InsetChart.json")).read())
        adult_vectors[iseed] = ic_json["Channels"]["Adult Vectors"]["Data"][-365:]

    adult_vector_means[isim] = adult_vectors.mean(axis=1)

##
# ------ Set fig paths
fig_dir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\presentations_writeups\\gene_drive_paper\\figures'
os.makedirs(fig_dir, exist_ok=True)

##
# ---- Plot adult vectors, infectious vectors (frac), infectious vectors (number), daily EIR
colors = [(87/255, 69/255, 93/255), (163/255, 33/255, 109/255), (97/255, 130/255, 193/255)]
alpha = 0.8

fig, ax = plt.subplots(1, 1, figsize=(12, 6))

# - Filled time series
for ieir in range(0, len(adult_vector_means[::-1])):
    if ieir == 0:
        ax.fill_between(list(range(0, 365)), adult_vector_means[ieir],
                        label=eir_labels[ieir], color=colors[ieir], alpha=alpha)
    elif ieir > 0:
        ax.fill_between(list(range(0, 365)), adult_vector_means[ieir-1], adult_vector_means[ieir],
                        label=eir_labels[ieir], color=colors[ieir], alpha=alpha)

# - Non-filled time series
# for ieir, adult_vector_mean in enumerate(adult_vector_means[::-1]):
#     ax.plot(list(range(0, 365)), adult_vector_mean,
#             label=eir_labels[ieir], lw=3)

# - Release day vertical line
# ax.axvline(x=released_day, color='k', linestyle='--')

ax.set_ylim([0, 8900])
ax.set_xlim([0, 365])
ax.set_ylabel('Adult Vectors')
ax.yaxis.set_ticklabels([])
ax.set_xticks([0.0, 30.417, 60.833, 91.25, 121.667, 152.083,
               182.5, 212.917, 243.333, 273.75, 304.167, 334.583])
ax.set_xticklabels(['Jan 1', 'Feb 1', 'Mar 1', 'Apr 1', 'May 1', 'Jun 1',
                    'Jul 1', 'Aug 1', 'Sep 1', 'Oct 1', 'Nov 1', 'Dec 1'],
                   rotation=45)
# ax.legend(loc='upper left', frameon=False)
# - Or reverse order of handles/labels if desired
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='upper left', frameon=False)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

fig.tight_layout()
fig_file_png = os.path.join(fig_dir, 'vector_numbers_seasonality.png')
fig_file_pdf = os.path.join(fig_dir, 'vector_numbers_seasonality.pdf')
plt.savefig(fig_file_png, bbox_inches="tight", dpi=300)
plt.savefig(fig_file_pdf, bbox_inches="tight", dpi=300)
plt.show()
