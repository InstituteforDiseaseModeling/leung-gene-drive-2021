import configparser
import gepi.grid.construction as grid
from input_file_generation.DemographicsGenerator import DemographicsGenerator
from input_file_generation.MigrationGenerator import MigrationGenerator
import json
import logging
import numpy as np
import os
import pandas as pd
from simtools.SetupParser import SetupParser


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
    code_block = 'LOCAL'
    SetupParser(code_block)

    country = 'Burkina_Faso'
    point_records_path_main = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\' \
                              'DTK_input_output_staging\\input\\BurkinaFasoSpatialInsideNodes'
    if not os.path.exists(point_records_path_main):
        os.mkdir(point_records_path_main)

    for i in range(0, 1):
        point_records_path = point_records_path_main
        if not os.path.exists(point_records_path):
            os.mkdir(point_records_path)
        create_grid_file = 0  # must be 1 when create_demographics_file == 1
        create_demographics_file = 0
        create_human_migration_file = 0
        create_vector_migration_file = 1

        if create_grid_file:
            # - Create grid from household lat-lons
            pop_filename = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\DTK_input_output_staging\\' \
                           'input\\population_data\\facebook_pop_clipped_BurkinaFasoSpatialInsideNodesAnd7kmBuffer.csv'
            if not os.path.exists(os.path.join(point_records_path, 'Grid', 'full_pop_grid.csv')):
                    final_grid_file = create_grid_files(point_records_path, pop_filename, 'full_pop')
            else:
                final_grid_file = os.path.join(point_records_path, 'Grid', 'full_pop_grid.csv')

        if create_demographics_file:
            # - Demographics filepath
            demo_files_dir = os.path.join(point_records_path, 'Demographics')
            if not os.path.exists(demo_files_dir):
                os.mkdir(demo_files_dir)

            demo_fp = os.path.join(demo_files_dir, "demographics.json")

            # - Generate demographics file
            demo = DemographicsGenerator.from_grid_file(population_input_file=final_grid_file,
                                                        demographics_filename=demo_fp)

        if create_human_migration_file:
            # - Demographics filepath
            demo_files_dir = os.path.join(point_records_path, 'Demographics')
            demo_fp = os.path.join(demo_files_dir, "demographics.json")
            if not os.path.exists(demo_fp):
                raise ValueError('Please provide a valid demographics file.')

            # - Migration directory
            human_migration_files_dir = os.path.join(point_records_path, 'HumanMigration', 'gm_def')
            if not os.path.exists(human_migration_files_dir):
                os.mkdir(human_migration_files_dir)
            mig_types = ['local']

            # - Generate human migration files
            for mig_type in mig_types:
                migration_fp = "human_%s_migration.bin" % (mig_type)

                # --> use original MigrationGenerator for this
                MigrationGenerator(demographics_file_path=demo_fp,
                                   migration_dir=human_migration_files_dir,
                                   migration_filename=migration_fp,
                                   gravity_migr_params=np.array(
                                       [7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00]),
                                   site="human",
                                   mig_type=mig_type).generate_migration(outfilename=migration_fp + '.json')

        if create_vector_migration_file:
            # - Demographics filepath
            demo_files_dir = os.path.join(point_records_path, 'Demographics')
            demo_fp = os.path.join(demo_files_dir, "demographics.json")
            if not os.path.exists(demo_fp):
                raise ValueError('Please provide a valid demographics file.')

            # - Migration directory
            vector_migration_files_dir = os.path.join(point_records_path, 'VectorMigration', 'thomas2013negexpovnn_b05')
            if not os.path.exists(vector_migration_files_dir):
                os.mkdir(vector_migration_files_dir)
            mig_types = ['local']

            # - Generate vector migration files
            for mig_type in mig_types:
                migration_fp = "vector_%s_migration.bin" % (mig_type)

                # --> use modified MigrationGenerator for this (modified_migration_generator.py)
                MigrationGenerator(demographics_file_path=demo_fp,
                                   migration_dir=vector_migration_files_dir,
                                   migration_filename=migration_fp,
                                   gravity_migr_params=np.array(
                                       [7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00]),  # not used
                                   site="vector",
                                   mig_type=mig_type).generate_migration(outfilename=migration_fp + '.json')
