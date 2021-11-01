import csv
import json
import pandas as pd
import numpy as np
from geopy.distance import geodesic
from struct import pack
import os
import dtk.tools.demographics.compiledemog as compiledemog
import re

from input_file_generation.create_migration_header import main

from simtools.Utilities.General import init_logging
logger = init_logging('MigrationGenerator')


class MigrationGenerator:
    """
    Generates migration file based on demographics file.
    """

    def __init__(self, demographics_file_path, migration_dir, migration_filename, gravity_migr_params, site, mig_type, **kwargs):

        self.demographics_file_path = demographics_file_path
        self.migration_dir = migration_dir
        self.migration_filename = migration_filename
        self.gravity_migr_params = gravity_migr_params
        self.site = site
        self.mig_type=mig_type
        self.kwargs = kwargs

    def load_demo(self, demo_file):
        with open(demo_file, 'r') as f:
            demo_dict = json.load(f)

        N = demo_dict['Metadata']['NodeCount']
        lat = np.ones(N)
        long = np.ones(N)
        grid_id = np.ones(N)
        node_id = np.ones(N)
        pop = np.ones(N)

        for i in range(N):
            node = demo_dict['Nodes'][i]
            lat[i] = node['NodeAttributes']['Latitude']
            long[i] = node['NodeAttributes']['Longitude']
            grid_id[i] = node['NodeAttributes']['FacilityName']
            node_id[i] = node['NodeID']
            pop[i] = node['NodeAttributes']['InitialPopulation']

        df = pd.DataFrame({
            'lat': lat,
            'long': long,
            'grid_id': grid_id,
            'node_id': node_id,
            'pop': pop
        })

        return df

    def compute_migr_prob(self, grav_params, ph, pd, d):
        # If home/dest node has 0 pop, assume this node is the regional work node-- no local migration allowed
        if ph == 0 or pd == 0:
            return 0.
        else:
            num_trips = grav_params[0] * ph ** grav_params[1] * pd ** grav_params[2] * d ** grav_params[3]
            prob_trip = np.min([1., num_trips / ph])
            return prob_trip

    def compute_migr_dict(self, df, grav_params, return_prob_sums=False, **kwargs):
        '''
        save link rates to a human readable file;
        the txt file is consumable by link_rates_txt_2_bin(self) like function to generate DTK migration binary
        '''
        migr = {}

        # ----------- ADDED
        # - Create array of discrete distances
        max_grid_dim = 20  # km; 20 km wide, 14 km high
        discrete_dists = np.array([])
        discrete_max_dyxs = np.array([])
        discrete_num_nodess = np.array([])
        # the following assumes 1 km grid resolution
        for idist1 in range(1, max_grid_dim+1):
            discrete_dists = np.append(discrete_dists, idist1)
            discrete_max_dyxs = np.append(discrete_max_dyxs, idist1)
            discrete_num_nodess = np.append(discrete_num_nodess, 4)
            for idist2 in range(1, idist1+1):
                discrete_dists = np.append(discrete_dists, np.sqrt(idist2**2 + idist1**2))
                discrete_max_dyxs = np.append(discrete_max_dyxs, idist1)
                if idist2 == idist1:
                    discrete_num_nodess = np.append(discrete_num_nodess, 4)
                else:
                    discrete_num_nodess = np.append(discrete_num_nodess, 8)
        # ----------- /ADDED

        p_sum = np.zeros(len(df))
        jj = 0
        for i1, r1 in df.iterrows():
            migr[r1['node_id']] = {}

            for i2, r2 in df.iterrows():
                if r2['node_id'] == r1['node_id']:
                    pass
                else:
                    d = geodesic((r1['lat'], r1['long']), (r2['lat'], r2['long'])).km

                    # ----------- ADDED
                    dy = geodesic((r1['lat'], r1['long']), (r2['lat'], r1['long'])).km
                    dx = geodesic((r1['lat'], r1['long']), (r1['lat'], r2['long'])).km
                    discrete_dist = find_nearest_discrete_dist(d, discrete_dists)
                    discrete_max_dyx = np.rint(np.max([dy, dx]))
                    discrete_num_nodes = np.asscalar(discrete_num_nodess[(discrete_dists == discrete_dist)
                                                                         & (discrete_max_dyxs == discrete_max_dyx)])

                    # - Thomas et al. (2013) after his corrections (see email)
                    # --neg exp model
                    b = 0.5
                    migr[r1['node_id']][r2['node_id']] = np.exp(-b * d) / discrete_num_nodes
                    # ----------- /ADDED

                    # - original gravity model
                    # migr[r1['node_id']][r2['node_id']] = self.compute_migr_prob(grav_params, r1['pop'], r2['pop'], d)

            p_sum[jj] = np.sum(list(migr[r1['node_id']].values()))
            print(p_sum[jj])
            jj += 1

        if return_prob_sums:
            return [migr, p_sum]
        else:
            return migr

    def gen_gravity_links_json(self, demo_file, grav_params, outf=None, **kwargs):
        df = self.load_demo(demo_file)

        if 'exclude_nodes' in kwargs:
            df = df[np.logical_not(np.in1d(df['node_id'],kwargs['exclude_nodes']))]

        migr_dict = self.compute_migr_dict(df, grav_params, return_prob_sums=False, **kwargs)

        # Save to file:
        if outf == None:
            outf = 'grav_migr_rates.json'
        with open(outf, 'w') as f:
            json.dump(migr_dict, f, indent=4)

        return migr_dict

    @staticmethod
    def save_link_rates_to_txt(rates_txt_file_path, link_rates):
        '''
        convert a txt links rates file (e.g. as generated by save_link_rates_to_txt(self...)) to DTK binary migration file
        '''
        with open(rates_txt_file_path, 'w') as fout:
            for src, v in link_rates.items():
                for dest, mig in v.items():
                    fout.write('%d %d %0.1g\n' % (int(src), int(dest), mig))

    def link_rates_txt_2_bin(self, rates_txt_file_path, rates_bin_file_path, route="local"):

        fopen = open(rates_txt_file_path)
        fout = open(rates_bin_file_path, 'wb')

        net = {}
        net_rate = {}

        MAX_DESTINATIONS_BY_ROUTE = {'local': 100,
                                     'regional': 30,
                                     'sea': 5,
                                     'air': 60}

        for line in fopen:
            s = line.strip().split()
            ID1 = int(float(s[0]))
            ID2 = int(float(s[1]))
            rate = float(s[2])
            # print(ID1,ID2,rate)
            if ID1 not in net:
                net[ID1] = []
                net_rate[ID1] = []
            net[ID1].append(ID2)
            net_rate[ID1].append(rate)

        for ID in sorted(net.keys()):

            ID_write = []
            ID_rate_write = []

            if len(net[ID]) > MAX_DESTINATIONS_BY_ROUTE[route]:
                print('There are %d destinations from ID=%d.  Trimming to %d (%s migration max) with largest rates.' % (
                    len(net[ID]), ID, MAX_DESTINATIONS_BY_ROUTE[route], route))
                dest_rates = list(zip(net[ID], net_rate[ID]))
                dest_rates.sort(key=lambda tup: tup[1], reverse=True)
                trimmed_rates = dest_rates[:MAX_DESTINATIONS_BY_ROUTE[route]]
                # print(len(trimmed_rates))
                # print(f"{trimmed_rates}")
                (net[ID], net_rate[ID]) = zip(*trimmed_rates)
                # print(net[ID], net_rate[ID])

            for i in range(MAX_DESTINATIONS_BY_ROUTE[route]):
                ID_write.append(0)
                ID_rate_write.append(0)
            for i in range(len(net[ID])):
                ID_write[i] = net[ID][i]
                ID_rate_write[i] = net_rate[ID][i]
            # print(ID_write)
            # print(ID_rate_write)
            s_write = pack('L' * len(ID_write), *ID_write)
            s_rate_write = pack('d' * len(ID_rate_write), *ID_rate_write)
            fout.write(s_write)
            fout.write(s_rate_write)

        fopen.close()
        fout.close()

    def save_migration_header(self, outfilename=None):

        # generate migration header for DTK consumption
        # todo: the script below needs to be refactored/rewritten
        # in its current form it requires compiled demographisc file (that's not the only problem with its design)
        # to compile the demographics file need to know about compiledemog file here, which is unnecessary
        # compiledemog.py too could be refactored towards object-orientedness
        # the demographics_file_path supplied here may be different from self.demographics_file_path)
        # mig_type = 'Local', 'Regional'
        compiledemog.main(self.demographics_file_path)
        main('dtk-tools', re.sub('\.json$', '.compiled.json', self.demographics_file_path), self.mig_type, outfilename=outfilename)

    def generate_migration(self, outfilename=None):
        """
        save migration binaries and migration header
        """
        migr_json_fp = os.path.join(self.migration_dir, "%s_grav_migr_%s_rates.json" % (self.site, self.mig_type))

        if 'exclude_nodes_from_regional_migration' in self.kwargs:
            exclude_nodes = self.kwargs['exclude_nodes_from_regional_migration']
            migr_dict = self.gen_gravity_links_json(self.demographics_file_path, self.gravity_migr_params, outf=migr_json_fp,
                                                    exclude_nodes=exclude_nodes)
        else:
            migr_dict = self.gen_gravity_links_json(self.demographics_file_path, self.gravity_migr_params, outf=migr_json_fp)

        rates_txt_fp = os.path.join(self.migration_dir, "%s_grav_migr_%s_rates.txt" % (self.site, self.mig_type))

        self.save_link_rates_to_txt(rates_txt_fp, migr_dict)

        # Generate migration binary:
        print("migration_filename: ", self.migration_filename)
        self.link_rates_txt_2_bin(rates_txt_fp, os.path.join(self.migration_dir, self.migration_filename))

        self.save_migration_header(outfilename=os.path.join(self.migration_dir, outfilename))


# ----------- ADDED
def find_nearest_discrete_dist(cont_dist, discrete_dists):
    idx = (np.abs(discrete_dists - cont_dist)).argmin()
    return discrete_dists[idx]
# ----------- /ADDED
