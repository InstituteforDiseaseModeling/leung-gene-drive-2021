import json
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rcParams
import os
import pandas as pd
import shapefile

rcParams['pdf.fonttype'] = 42

# ------ Set fig path
fig_dir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\presentations_writeups\\gene_drive_paper\\figures'
os.makedirs(fig_dir, exist_ok=True)

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
# ------ Load spatial plot shapefiles + set map boundaries

# - Set map boundaries
lngmin = min(dfg['lon']) - 0.02
lngmax = max(dfg['lon']) + 0.02
latmin = min(dfg['lat']) - 0.02
latmax = max(dfg['lat']) + 0.02
minmaxes = (lngmin, lngmax, latmin, latmax)


# - Load roads and rivers
def features_within_polygon(shapefile_path, lngmin, lngmax, latmin, latmax):
    """
    Plot man-made features like roads or natural features like rivers
    :param shapefile_path: path to shapefile describing feature
    :param lngmin: westernmost longitude of polygon to plot
    :param lngmax: easternmost longitude of polygon to plot
    :param latmin: southernmost longitude of polygon to plot
    :param latmax: northernmost longitude of polygon to plot
    :return: features: coordinates of features
    """

    sf = shapefile.Reader(shapefile_path)
    features = []
    for shape in sf.shapeRecords():
        x = [i[0] for i in shape.shape.points[:]]
        x2 = [i for i in x if lngmin < i < lngmax]
        y = [i[1] for i in shape.shape.points[:]]
        y2 = [i for i in y if latmin < i < latmax]
        if x2 and y2:
            features.append(shape.shape.points[:])

    return features


river_shapefile_path = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)',
                                    'Malaria Team Folder', 'data',
                                    'Burkina', 'Burkina shapefiles', 'BFA_wat', 'BFA_water_lines_dcw')
road_shapefile_path = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)',
                                   'Malaria Team Folder', 'data',
                                   'Burkina', 'Burkina shapefiles', 'BFA_rds', 'BFA_roads')
rivers = features_within_polygon(river_shapefile_path, lngmin=lngmin, lngmax=lngmax, latmax=latmax, latmin=latmin)
roads = features_within_polygon(road_shapefile_path, lngmin=lngmin, lngmax=lngmax, latmax=latmax, latmin=latmin)

##
# ------ Plot population on spatial grid

# - Set color map
cmap = cm.viridis

fig, ax = plt.subplots(1, 1, figsize=(12, 12))

# - Plot rivers
for river in rivers:
    x = [i[0] for i in river]
    y = [i[1] for i in river]
    ax.plot(x, y, zorder=1, color='#7faddd', lw=5)

# - Plot roads
for road in roads:
    x = [i[0] for i in road]
    y = [i[1] for i in road]
    ax.plot(x, y, zorder=1, color='#969696', lw=3)

ax.set_ylim([minmaxes[2], minmaxes[3]])
ax.set_xlim([minmaxes[0], minmaxes[1]])

# - Color grid cells by population
scatter_size = 690
sc = ax.scatter(dfg['lon'], dfg['lat'], zorder=2, marker='s',
                s=scatter_size, c=dfg['pop'], cmap=cmap, edgecolor='silver')
cb = fig.colorbar(sc, ax=[ax], fraction=0.033, pad=0.025, location='left')
cb.set_label(label='Human population', size=18)
for t in cb.ax.get_yticklabels():
     t.set_fontsize(18)

# - Highlight gene drive release nodes
dfg_rn = dfg.sort_values(by=['pop'], ascending=False)
dfg_rn = dfg_rn[0:6]
ax.scatter(dfg_rn['lon'], dfg_rn['lat'], zorder=2, marker='s',
           s=scatter_size-20, facecolor='None', edgecolor='red', linewidth=3)

# - Annotate each grid cell w/ pop numbers
# for irow in range(0, len(dfg)):
#     dfnow = dfg.iloc[irow]
#     ax.text(dfnow['lon'], dfnow['lat'], str(dfnow['pop']), fontsize=10, fontweight='bold', color='silver')

# - Plot scalebar
ax.plot([-1.57267, -1.52673], [12.02, 12.02], 'k', lw=2)
ax.text(-1.57267, 12.012, r'5 km', fontsize=21)

# - Make plot square
ax.set_aspect('equal', adjustable='box')

# - Remove x and y axis labels
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)

# - Save plot
fig_file_png = os.path.join(fig_dir, 'spatialinside_human_pop.png')
fig_file_pdf = os.path.join(fig_dir, 'spatialinside_human_pop.pdf')
plt.savefig(fig_file_png, bbox_inches="tight", dpi=300)
plt.savefig(fig_file_pdf, bbox_inches="tight", dpi=300)
plt.show()

##
# ------ Calculate total and release nodes populations
print(dfg['pop'].sum())
# 3670
print(dfg_rn['pop'].sum())
# 833
