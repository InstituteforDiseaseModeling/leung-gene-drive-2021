# ----------------------------------------------------------------------
# Import and set universal params
# ----------------------------------------------------------------------
##
from geopy.distance import geodesic
import pandas as pd
import numpy as np
import os
import json
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42

# ----------------------------------------------------------------------
# Set paths
# ----------------------------------------------------------------------
##
# ------ Set data paths
# - spatialinside_classic3allele_no_intvns_aEIR30_vecmig_sweep_rr0, rr0 = 0.1,
# migration reported: start_day=212, end_day=274  # aug 1 through oct 1
migration_mul, rn, vmig_dir = 1, 0, 'thomas2013negexpovnn_b05'
data_dir = 'Y:\\home\\sleung\\output\\spatialinside_classic3allele_no_in_20211013_064338\\5ca\\4a9\\eff\\5ca4a9ef-f02b-ec11-9ecd-9440c9bee941'

# ------ Set fig paths
fig_dir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\presentations_writeups\\gene_drive_paper\\figures'
os.makedirs(fig_dir, exist_ok=True)

# ----------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------
##
# ------ Load vector data
filem = os.path.join(data_dir, 'output\\ReportVectorMigration.csv')
dfm = pd.read_csv(filem)
if 'Unnamed: 0' in dfm.columns:
    dfm = dfm.drop('Unnamed: 0', axis=1)

files = os.path.join(data_dir, 'output\\ReportVectorStats.csv')
dfs = pd.read_csv(files)
if 'Unnamed: 0' in dfs.columns:
    dfs = dfs.drop('Unnamed: 0', axis=1)

##
# ------ Load grid/demographics data
demographics_dir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\' \
                   'DTK_input_output_staging\\input\\BurkinaFasoSpatialInsideNodes\\Demographics\\'
demographics_file = 'demographics.json'
demographics_filepath = os.path.join(demographics_dir, demographics_file)
nodes = []
pop = []
lat = []
lon = []
with open(demographics_filepath) as f:
    data = json.load(f)
    for node in data['Nodes']:
        nodes.append(node['NodeID'])
        pop.append(node['NodeAttributes']['InitialPopulation'])
        lat.append(node['NodeAttributes']['Latitude'])
        lon.append(node['NodeAttributes']['Longitude'])
dfg = pd.DataFrame({'NodeID': nodes, 'pop': pop, 'lat': lat, 'lon': lon})

##
# ------ Clean data
dfm.drop([' MigrationType', ' Species', ' Age'], axis=1, inplace=True)
dfm.rename(columns={' ID': 'ID',
                    ' FromNodeID': 'FromNodeID',
                    ' ToNodeID': 'ToNodeID'}, inplace=True)

# ----------------------------------------------------------------------
# Calc distance btwn nodes
# ----------------------------------------------------------------------
##
# ------ Create array of discrete distances
# HARDCODE
max_grid_dim = 20  # km; 20 km wide, 14 km high
discrete_dists = np.array([0])
discrete_max_dyxs = np.array([0])
discrete_num_nodess = np.array([0])
# the following assumes 1 km grid resolution
for idist1 in range(1, max_grid_dim + 1):
    discrete_dists = np.append(discrete_dists, idist1)
    discrete_max_dyxs = np.append(discrete_max_dyxs, idist1)
    discrete_num_nodess = np.append(discrete_num_nodess, 4)
    # print('idist1: ' + str(idist1))
    for idist2 in range(1, idist1 + 1):
        discrete_dists = np.append(discrete_dists, np.sqrt(idist2 ** 2 + idist1 ** 2))
        # print('idist2: ' + str(idist2))
        discrete_max_dyxs = np.append(discrete_max_dyxs, idist1)
        if idist2 == idist1:
            discrete_num_nodess = np.append(discrete_num_nodess, 4)
        else:
            discrete_num_nodess = np.append(discrete_num_nodess, 8)


##
# ------ Define fxns
def distance_from_to_nodeID(fromNodeID, toNodeID):
    from_lat = dfg[dfg['NodeID'] == fromNodeID]['lat'].iloc[0]
    from_lon = dfg[dfg['NodeID'] == fromNodeID]['lon'].iloc[0]
    to_lat = dfg[dfg['NodeID'] == toNodeID]['lat'].iloc[0]
    to_lon = dfg[dfg['NodeID'] == toNodeID]['lon'].iloc[0]
    d_km = geodesic((from_lat, from_lon), (to_lat, to_lon)).km
    return d_km


def find_nearest_discrete_dist(row):
    idx = (np.abs(discrete_dists - row['Distance'])).argmin()
    return discrete_dists[idx]


def find_discrete_max_dyx(fromNodeID, toNodeID):
    from_lat = dfg[dfg['NodeID'] == fromNodeID]['lat'].iloc[0]
    from_lon = dfg[dfg['NodeID'] == fromNodeID]['lon'].iloc[0]
    to_lat = dfg[dfg['NodeID'] == toNodeID]['lat'].iloc[0]
    to_lon = dfg[dfg['NodeID'] == toNodeID]['lon'].iloc[0]
    dy = geodesic((from_lat, from_lon), (to_lat, from_lon)).km
    dx = geodesic((from_lat, from_lon), (from_lat, to_lon)).km
    discrete_max_dyx = np.rint(np.max([dy, dx]))
    return discrete_max_dyx


##
# ------ Set node characteristics
calcnow = 1
if calcnow == 1:

    dbndf = pd.DataFrame(columns=['FromNodeID', 'ToNodeID', 'Distance', 'Max Distxy'])

    for fromnode in nodes:
        for tonode in nodes:
            dbndf = dbndf.append({'FromNodeID': fromnode, 'ToNodeID': tonode,
                                  'Distance': distance_from_to_nodeID(fromnode, tonode),
                                  'Max Distxy': find_discrete_max_dyx(fromnode, tonode)},
                                 ignore_index=True)

##
# ------ Calc distance btwn all migrations
dfm = pd.merge(dfm, dbndf, how='left', on=['FromNodeID', 'ToNodeID'])
dfm['Discrete Distance'] = dfm.apply(find_nearest_discrete_dist, axis=1)

# ----------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------
##
# ------ Plot total # of migrations vs. distance (over reported number of days)
# (#1 plot type in notebook for 2/2/21)
plotnow = 1
if plotnow == 1:
    max_bin_edge = np.ceil(dfm['Distance'].max()) + 1

    fig, axes = plt.subplots(2, 1, figsize=(12, 12))
    for iax, ax in enumerate(axes):
        if iax == 0:
            ax.hist(dfm['Distance'], bins=np.linspace(0, max_bin_edge, 2 * max_bin_edge + 1))
            ax.set_ylabel('# of migrations')
        elif iax == 1:
            hist, bins = np.histogram(dfm['Distance'], bins=np.linspace(0, max_bin_edge, 2 * max_bin_edge + 1))
            ax.bar(bins[:-1] + 0.25, 100 * hist.astype(np.float32) / hist.sum(),
                   width=(bins[1] - bins[0]))
            ax.set_ylabel('% of migrations')
        ax.set_xlabel('Distance (km)')
        ax.set_title('migration_mul = ' + str(migration_mul))
    fig_file = os.path.join(fig_dir, 'hist0.5km_total_mig_vs_dist_rn' + str(rn)
                            + '_' + vmig_dir + '_vmigmul' + str(migration_mul) + '.png')
    migration_mul, rn, vmig_dir = 1, 0, 'thomas2013negexpovnn_b05'
    plt.savefig(fig_file, dpi=300)
    plt.show()

    fig, axes = plt.subplots(2, 1, figsize=(12, 12))
    for iax, ax in enumerate(axes):
        if iax == 0:
            ax.hist(dfm['Distance'], bins=np.linspace(0, max_bin_edge, max_bin_edge + 1))
            ax.set_ylabel('# of migrations')
        elif iax == 1:
            hist, bins = np.histogram(dfm['Distance'], bins=np.linspace(0, max_bin_edge, max_bin_edge + 1))
            ax.bar(bins[:-1] + 0.5, 100 * hist.astype(np.float32) / hist.sum(),
                   width=(bins[1] - bins[0]))
            ax.set_ylabel('% of migrations')
        ax.set_xlabel('Distance (km)')
        ax.set_title('migration_mul = ' + str(migration_mul))
    fig_file = os.path.join(fig_dir, 'hist1km_total_mig_vs_dist_rn' + str(rn)
                            + '_' + vmig_dir + '_vmigmul' + str(migration_mul) + '.png')
    plt.savefig(fig_file, dpi=300)
    plt.show()

##
# ------ Calculate individual mosquito movements
calcnow = 1
if calcnow == 1:
    indiv_dist_means = np.full(len(dfm['ID'].unique()), np.nan)
    indiv_dist_stds = np.full(len(dfm['ID'].unique()), np.nan)
    indiv_dist_sums = np.full(len(dfm['ID'].unique()), np.nan)
    indiv_mvs = np.full(len(dfm['ID'].unique()), np.nan)

    for ivid, vid in enumerate(dfm['ID'].unique()):
        dfmnow = dfm[dfm['ID'] == vid]
        indiv_dist_means[ivid] = dfmnow['Distance'].mean()
        indiv_dist_stds[ivid] = dfmnow['Distance'].std()
        indiv_dist_sums[ivid] = dfmnow['Distance'].sum()
        indiv_mvs[ivid] = len(dfmnow)

##
# ------ Plot mean and std of individual mosq movements
plotnow = 1
if plotnow == 1:
    fig, axes = plt.subplots(5, 1, figsize=(11, 25))

    ax = axes[0]
    ax.scatter(np.arange(0, len(indiv_dist_means)), indiv_dist_means, s=1, alpha=0.5)
    ax.set_xlabel('Individual vector')
    ax.set_ylabel('Avg dist trvled in 1 ts [km]')

    ax = axes[1]
    ax.scatter(np.arange(0, len(indiv_dist_means)), indiv_dist_sums, s=1, alpha=0.5)
    ax.set_xlabel('Individual vector')
    ax.set_ylabel('Summed dist trvled [km]')

    ax = axes[2]
    hist, bins = np.histogram(indiv_dist_means, bins=np.linspace(0, max_bin_edge, 2 * max_bin_edge + 1))
    ax.bar(bins[:-1] + 0.25, hist, width=0.5)
    ax.set_xlabel('Avg dist trvled in 1 ts [km]')
    ax.set_ylabel('# of indiv vectors')

    ax = axes[3]
    hist, bins = np.histogram(indiv_dist_sums, bins=50)
    ax.bar(bins[:-1] + (bins[1] - bins[0])/2, hist,
           width=bins[1] - bins[0])
    ax.set_xlabel('Summed dist trvled [km]')
    ax.set_ylabel('# of indiv vectors')

    ax = axes[4]
    nummoves, numnummoves = np.unique(indiv_mvs, return_counts=True)
    ax.scatter(nummoves, numnummoves)
    for i, v in enumerate(numnummoves):
        ax.annotate(str(v), xy=(nummoves[i], v), xytext=(-7, 7), textcoords='offset points')
    ax.set_xlabel('# of moves')
    ax.set_ylabel('# of indiv vectors')

    fig_file = os.path.join(fig_dir, 'indiv_vector_distances_rn' + str(rn)
                            + '_' + vmig_dir + '_vmigmul' + str(migration_mul) + '.png')
    plt.savefig(fig_file, dpi=300)
    plt.show()
