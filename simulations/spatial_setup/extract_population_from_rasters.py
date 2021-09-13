import os
import rasterio
import numpy as np
import pandas as pd
from affine import Affine
from pyproj import Proj, transform
import shapefile
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

from raster_clipping_and_extraction_functions import run_all_clipping
from plotting.colors import load_color_palette

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams.update({'font.size': 20})


def extract_blob_coordinates(shapefile_path):
    sf = shapefile.Reader(shapefile_path)
    features = []
    sns.set_style('whitegrid', {'axes.linewidth': 0.5})
    for shape in sf.shapeRecords():
        features.append(shape.shape.points[:])

    return features


def features_within_polygon(shapefile_path, lngmin, lngmax, latmin, latmax):
    '''
    Plot man-made features like roads or natural features like rivers
    :param shapefile_path: path to shapefile describing feature
    :param lngmin: westernmost longitude of polygon to plot
    :param lngmax: easternmost longitude of polygon to plot
    :param latmin: southernmost longitude of polygon to plot
    :param latmax: northernmost longitude of polygon to plot
    :return: features: coordinates of features
    '''

    sf = shapefile.Reader(shapefile_path)
    features = []
    sns.set_style('whitegrid', {'axes.linewidth': 0.5})
    for shape in sf.shapeRecords():
        x = [i[0] for i in shape.shape.points[:]]
        x2 = [i for i in x if lngmin < i < lngmax]
        y = [i[1] for i in shape.shape.points[:]]
        y2 = [i for i in y if latmin < i < latmax]
        if x2 and y2:
            features.append(shape.shape.points[:])

    return features


def convert_raster_to_csv(fname):
    with rasterio.open(fname) as r:
        T0 = r.transform  # upper-left pixel corner affine transform
        p1 = Proj(r.crs)
        A = r.read()  # pixel values

    # All rows and columns
    cols, rows = np.meshgrid(np.arange(A.shape[2]), np.arange(A.shape[1]))

    # Get affine transform for pixel centres
    T1 = T0 * Affine.translation(0.5, 0.5)
    # Function to convert pixel row/column index (from 0) to easting/northing at centre
    rc2en = lambda r, c: (c, r) * T1

    # All eastings and northings (there is probably a faster way to do this)
    eastings, northings = np.vectorize(rc2en, otypes=[np.float, np.float])(rows, cols)

    # Project all longitudes, latitudes
    p2 = Proj(proj='latlong', datum='WGS84')
    longs, lats = transform(p1, p2, eastings, northings)
    values = A[0]
    values[values < 0] = 0
    values = list(values.ravel())
    lats = list(lats.ravel())
    longs = list(longs.ravel())
    df = pd.DataFrame({'latitude': lats, 'longitude': longs, 'pop': values})
    df.to_csv(fname.split('.')[0] + '.csv')

    return df


def plot_data_with_features(data, rivers, roads, blob_coords, minmaxes, figpath, **kwargs):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.gca()

    blob_colors = ['r', 'g', 'k']
    for j, blobs in enumerate(blob_coords):
        for blob in blobs:
            x = [i[0] for i in blob]
            y = [i[1] for i in blob]
            ax.plot(x, y, color=blob_colors[j])

    for river in rivers:
        x = [i[0] for i in river]
        y = [i[1] for i in river]
        ax.plot(x, y, color='#7faddd')

    for road in roads:
        x = [i[0] for i in road]
        y = [i[1] for i in road]
        ax.plot(x, y, color='#969696')

    data = data[data['pop'] != 0]
    ax.scatter(data['longitude'], data['latitude'], 20, linewidth=0,
               color='#969696',
               alpha=0.5,
               label='not selected')
    ax.set_ylim([minmaxes[2], minmaxes[3]])
    ax.set_xlim([minmaxes[0], minmaxes[1]])

    if 'indie' in kwargs:
        palette = load_color_palette()
        palette = [palette[x] for x in [4, 5, 9]]

        for i, arm in enumerate(['arm1', 'arm2', 'arm3']):
            column = 'pop_selected_' + arm
            data_selected = data[data[column] != 0]
            ax.scatter(data_selected['longitude'], data_selected['latitude'], 20, linewidth=0,
                       color=palette[i],
                       label=arm)
        ax.legend()

    plt.savefig(figpath)

    return None


if __name__ == '__main__':

    # set blob to any whole number not between 1 and 3 to plot entire area

    blob = 0

    # Paths to roads and rivers shape file, and to full Burkina Faso raster file
    river_shapefile_path = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)',
                                        'Malaria Team Folder', 'data',
                                        'Burkina', 'Burkina shapefiles', 'BFA_wat', 'BFA_water_lines_dcw')

    road_shapefile_path = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)',
                                       'Malaria Team Folder', 'data',
                                       'Burkina', 'Burkina shapefiles', 'BFA_rds', 'BFA_roads')

    indie_data_path = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)',
                                   'Malaria Team Folder', 'data',
                                   'Burkina', 'INDIE', 'Total_HH_pop_GIS_Nov_2018.csv')

    # out_path = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)',
    #                         'Malaria Team Folder', 'projects',
    #                         'Vector_genetics', 'Data')
    out_path = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\DTK_input_output_staging\\input\\population_data'

    burkina_raster_file_full = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)',
                                            'Malaria Team Folder', 'data',
                                            'Burkina', 'Burkina raster files', 'hrsl_bfa.tif')

    # grid_data = pd.read_csv(indie_data_path)
    grid_data0 = pd.read_csv(indie_data_path)
    bufferdist = 7  # km
    # ~111 km in 1 deg lat
    minlat = grid_data0['latitude'].min() - bufferdist/111
    maxlat = grid_data0['latitude'].max() + bufferdist/111
    # ~109 km in 1 deg lon at 12N: http://www.csgnetwork.com/degreelenllavcalc.html
    minlon = grid_data0['longitude'].min() - bufferdist/109
    maxlon = grid_data0['longitude'].max() + bufferdist/109
    grid_data = pd.DataFrame(
        {'latitude': [minlat, maxlat],
         'longitude': [minlon, maxlon]
         }
    )

    outfilename = 'facebook_pop_clipped_BurkinaFasoSpatialInsideNodesAnd7kmBuffer.csv'
    shp_path = out_path
    shapefile_type = 'bbox'

    if not os.path.exists(os.path.join(out_path, outfilename)):
        # Clip raster to suit lat-longs in grid data
        rasterpath = run_all_clipping(burkina_raster_file_full, out_path, shapefile_type, outfilename,
                                      grid_data, shp_path, overwrite=False,
                                      alpha=0.001, out_name="raster", write_shp=True, crop=True)

        # Convert facebook clipped raster data to csv similar to grid_data
        facebook_data = convert_raster_to_csv(rasterpath)
    else:
        facebook_data = pd.read_csv(os.path.join(out_path, outfilename))

    # if blob < 1 or blob > 3:
    #     # Get coordinates of roads and rivers in plots
    #     lngmin = min(grid_data['longitude'])
    #     lngmax = max(grid_data['longitude'])
    #     latmin = min(grid_data['latitude'])
    #     latmax = max(grid_data['latitude'])
    #     minmaxes = (lngmin, lngmax, latmin, latmax)
    #     rivers = features_within_polygon(river_shapefile_path, lngmin=lngmin, lngmax=lngmax, latmax=latmax,
    #                                      latmin=latmin)
    #     roads = features_within_polygon(road_shapefile_path, lngmin=lngmin, lngmax=lngmax, latmax=latmax, latmin=latmin)

    #     print("Total INDIE population: %f\n" % (sum(grid_data['hh_size'])))

    #     df_facebook = pd.read_csv(os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)',
    #                                            'Malaria Team Folder', 'projects',
    #                                            'Vector_genetics', 'Data',
    #                                            'facebook_pop_clipped.csv'))
    #     print("Total facebook population: %f\n" % (sum(df_facebook['pop'])))
