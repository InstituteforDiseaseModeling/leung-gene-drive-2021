import os
import json
import logging
import random
import numpy as np

import pandas as pd

import configparser

import gepi.grid.construction as grid
from simtools.SetupParser import SetupParser

from input_file_generation.DemographicsGenerator import DemographicsGenerator
from input_file_generation.MigrationGenerator import MigrationGenerator
from input_file_generation.ClimateGenerator import ClimateGenerator
# from input_file_generation.add_properties_to_demographics import generate_demographics_properties


def create_grid_files(path, records_filename, site):
    prefix = "%s_grid" % site
    point_records_path = path
    final_grid_files_dir = os.path.join(point_records_path, 'Grid')
    if not os.path.exists(final_grid_files_dir):
        os.mkdir(final_grid_files_dir)
    point_records_file = records_filename

    logging.info("Reading data...")
    point_records = pd.read_csv(os.path.join(point_records_file), encoding="iso-8859-1")
    point_records.rename(columns={'longitude': 'lon', 'latitude': 'lat'}, inplace=True)

    if not 'pop' in point_records.columns:
        point_records['pop'] = [5.5] * len(point_records)

    if 'hh_size' in point_records.columns:
        point_records['pop'] = point_records['hh_size']

    # point_records = point_records[point_records['pop']>0]
    x_min, y_min, x_max, y_max = grid.get_bbox(point_records)
    point_records = point_records[
        (point_records.lon >= x_min) & (point_records.lon <= x_max) & (point_records.lat >= y_min) & (
                point_records.lat <= y_max)]
    gridd, grid_id_2_cell_id, origin, final = grid.construct(x_min, y_min, x_max, y_max)
    gridd.to_csv(os.path.join(final_grid_files_dir, "%s_grid.csv" % site))

    with open(os.path.join(final_grid_files_dir, "%s_grid_id_2_cell_id.json" % site), "w") as g_f:
        json.dump(grid_id_2_cell_id, g_f, indent=3)

    point_records[['gcid', 'gidx', 'gidy']] = point_records.apply(grid.point_2_grid_cell_id_lookup,
                                                                  args=(grid_id_2_cell_id, origin,), axis=1).apply(
        pd.Series)

    grid_pop = point_records.groupby(['gcid', 'gidx', 'gidy'])['pop'].apply(np.sum).reset_index()
    grid_pop['pop'] = grid_pop['pop'].apply(lambda x: round(x / 5))
    grid_final = pd.merge(gridd, grid_pop, on='gcid')
    grid_final['node_label'] = list(grid_final.index)
    grid_final = grid_final[grid_final['pop'] > 5]
    grid_final.to_csv(os.path.join(final_grid_files_dir, prefix + '.csv'))

    return os.path.join(final_grid_files_dir, prefix + '.csv')


if __name__ == '__main__':

    config = configparser.ConfigParser()
    config.read('simtools.ini')
    # code_block = 'HPC'
    code_block = 'LOCAL'
    SetupParser(code_block)

    # site = "vector_genetics"  # prefix for generated demographics and migration files
    country = 'Burkina_Faso'
    # point_records_path_main = os.path.join(config['HPC']['input_root'])
    point_records_path_main =\
        'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\DTK_input_output_staging\\input\\BurkinaFasoSpatialInsideNodes'
    if not os.path.exists(point_records_path_main):
        os.mkdir(point_records_path_main)

    for i in range(0, 1):
        # point_records_path = os.path.join(point_records_path_main, 'Config_%i' %i)
        point_records_path = point_records_path_main
        if not os.path.exists(point_records_path):
            os.mkdir(point_records_path)
        create_grid_file = 0  # must be 1 when create_demographics_file == 1
        create_demographics_file = 0
        create_human_migration_file = 0
        create_vector_migration_file = 1
        create_climate_file = 0

        if create_grid_file:
            # - Create grid from household lat-lons
            # pop_filename = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)',
            #                             'Malaria Team Folder', 'projects',
            #                             'Vector_genetics', 'Data', 'facebook_pop_clipped.csv')
            pop_filename = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\DTK_input_output_staging\\input\\population_data\\facebook_pop_clipped_BurkinaFasoSpatialInsideNodesAnd7kmBuffer.csv'
            # if not os.path.exists(os.path.join(point_records_path, 'Grid', site + '_full_pop_grid.csv')):
            if not os.path.exists(os.path.join(point_records_path, 'Grid', 'full_pop_grid.csv')):
                    final_grid_file = create_grid_files(point_records_path, pop_filename, 'full_pop')
            else:
                # final_grid_file = os.path.join(point_records_path, 'Grid', site + '_full_pop_grid.csv')
                final_grid_file = os.path.join(point_records_path, 'Grid', 'full_pop_grid.csv')

        if create_demographics_file:
            # - Demographics filepath
            demo_files_dir = os.path.join(point_records_path, 'Demographics')
            if not os.path.exists(demo_files_dir):
                os.mkdir(demo_files_dir)

            # demo_fp = os.path.join(demo_files_dir, "%s_demographics.json" % site)
            demo_fp = os.path.join(demo_files_dir, "demographics.json")

            # - Generate demographics file
            demo = DemographicsGenerator.from_grid_file(population_input_file=final_grid_file,
                                                        demographics_filename=demo_fp)

        # ----- BEWARE THAT HUMAN AND VECTOR MIGRATION USE SAME MigrationGenerator
        # FXNS, CHANGE MigrationGenerator CODE IN BETWEEN CREATING HUMAN AND VECTOR
        # MIGRATION FILES!

        if create_human_migration_file:
            # - Demographics filepath
            demo_files_dir = os.path.join(point_records_path, 'Demographics')
            # demo_fp = os.path.join(demo_files_dir, "%s_demographics.json" % site)
            demo_fp = os.path.join(demo_files_dir, "demographics.json")
            if not os.path.exists(demo_fp):
                raise ValueError('Please provide a valid demographics file.')

            # - Migration directory
            # human_migration_files_dir = os.path.join(point_records_path, 'HumanMigration')
            human_migration_files_dir = os.path.join(point_records_path, 'HumanMigration', 'gm_def')
            if not os.path.exists(human_migration_files_dir):
                os.mkdir(human_migration_files_dir)
            mig_types = ['local']

            # - Generate human migration files
            for mig_type in mig_types:
                # migration_fp = "%s_%s_migration.bin" % (site, mig_type)
                migration_fp = "human_%s_migration.bin" % (mig_type)

                MigrationGenerator(demographics_file_path=demo_fp,
                                   migration_dir=human_migration_files_dir,
                                   migration_filename=migration_fp,
                                   gravity_migr_params=np.array(
                                       # [7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -3]),
                                       # [7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -2.5]),
                                       [7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00]),
                                   site="human",
                                   mig_type=mig_type).generate_migration(outfilename=migration_fp + '.json')

        if create_vector_migration_file:
            # - Demographics filepath
            demo_files_dir = os.path.join(point_records_path, 'Demographics')
            # demo_fp = os.path.join(demo_files_dir, "%s_demographics.json" % site)
            demo_fp = os.path.join(demo_files_dir, "demographics.json")
            if not os.path.exists(demo_fp):
                raise ValueError('Please provide a valid demographics file.')

            # - Migration directory
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovnn_b35')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovnn_b3')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovnn_b25')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovnn_b2')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovnn_b15')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovnn_b1')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b03')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b04')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b06')
            vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b07')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b027')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b032')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b037')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b042')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b047')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b005')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b01')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b025')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b05')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b075')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b1')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnnle_b01')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnnle_b025')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnnle_b05')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnnle_b075')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnnle_b1')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnnle')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnnle_b2')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnnle_b3')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnnle_b4')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnnle_b5')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovdle')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovdle_b0.2')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovdle_b0.3')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovdle_b0.4')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovdle_b0.5')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovdle_b0.6')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovdle_b0.7')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013halfcauchyovdle_b0.8')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovdle')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovdle_b1')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovdle_b2')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovdle_b3')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovdle_b4')
            # vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovdle_b5')
            if not os.path.exists(vector_migration_files_dir):
                os.mkdir(vector_migration_files_dir)
            mig_types = ['local']

            # - Generate vector migration files
            for mig_type in mig_types:
                # migration_fp = "%s_%s_migration.bin" % (site, mig_type)
                migration_fp = "vector_%s_migration.bin" % (mig_type)

                MigrationGenerator(demographics_file_path=demo_fp,
                                   migration_dir=vector_migration_files_dir,
                                   migration_filename=migration_fp,
                                   gravity_migr_params=np.array(
                                       [7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00]),
                                   site="vector",
                                   mig_type=mig_type).generate_migration(outfilename=migration_fp + '.json')

        if create_climate_file:
            # - Demographics filepath
            demo_files_dir = os.path.join(point_records_path, 'Demographics')
            # demo_fp = os.path.join(demo_files_dir, "%s_demographics.json" % site)
            demo_fp = os.path.join(demo_files_dir, "demographics.json")
            if not os.path.exists(demo_fp):
                raise ValueError('Please provide a valid demographics file.')

            # - Climate directory
            climate_files_dir = os.path.join(point_records_path, 'Climate')
            if not os.path.exists(climate_files_dir):
                os.mkdir(climate_files_dir)
            climate_logs_dir = os.path.join(point_records_path, 'Climate', 'Logs')
            if not os.path.exists(climate_logs_dir):
                os.mkdir(climate_logs_dir)

            # - Generate climate files
            cg = ClimateGenerator(demographics_file_path=demo_fp,
                                  work_order_path=os.path.join(climate_logs_dir, 'climate_wo.json'),
                                  climate_files_output_path=climate_files_dir,
                                  start_year=str(2016),
                                  num_years=str(2),
                                  climate_project="IDM-%s" % country,
                                  resolution=str(0),
                                  project_root='v2017').generate_climate_files()
