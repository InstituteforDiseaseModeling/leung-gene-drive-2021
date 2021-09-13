import numpy as np
import os
import json
import pandas as pd
import configparser

from dtk.utils.Campaign.CampaignClass import WaningEffectExponential
from dtk.interventions.mosquito_release import add_mosquito_release
from dtk.interventions.novel_vector_control import add_ATSB
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.vector.species import set_species_trait_modifiers, set_species_drivers, set_species_genes

from createSimDirectoryMapBR import createSimDirectoryMap


# --------- Define functions for reading serialized files
def get_outpath_for_serialized_file(simmap, comps_platform='Calculon',
                                    genes_type="classic",
                                    rr0=0, rr10=0, rr20=0, seed=0):

    temp = simmap[simmap['Run_Number'] == seed]

    if genes_type == "classic":
        temp = temp[temp['rr0'] == rr0]
    elif genes_type == "integral":
        temp = temp[(temp['rr10'] == rr10) & (temp['rr20'] == rr20)]

    if comps_platform == 'Calculon':
        outpath = temp['outpath'].values[0]
        outpath = outpath.replace('internal.idm.ctr', 'mnt')
        outpath = outpath.replace('\\', '/')
        outpath = outpath.replace('IDM2', 'idm2')
        return outpath
    elif comps_platform == 'Belegost':
        return temp['outpath'].values[0]


def make_simmap(expname):

    simmap = createSimDirectoryMap(expname)

    return simmap


# --------- Define fxns for vector migration
def update_vmig_mul(cb, vmig_mul):

    cb.update_params({
        'x_Vector_Migration_Local': vmig_mul  # vector migration
    })

    return {'vmig_mul': vmig_mul}


def update_vmig_dir(cb, input_file_dir, vmig_dir):

    if vmig_dir is not None:
        cb.update_params({
            'Vector_Migration_Filename_Local': os.path.join(input_file_dir, "VectorMigration", vmig_dir,
                                                            "vector_local_migration.bin"),
            'Enable_Vector_Migration': 1,
            'Enable_Vector_Migration_Local': 1,
            'Enable_Vector_Migration_Human': 0,
            'Enable_Vector_Migration_Wind': 0,
            'Vector_Sampling_Type': "TRACK_ALL_VECTORS",
            'Vector_Migration_Food_Modifier': 0,
            'Vector_Migration_Habitat_Modifier': 0,
            'Vector_Migration_Modifier_Equation': "LINEAR",
            'Vector_Migration_Stay_Put_Modifier': 0
        })

    return {'vmig_dir': vmig_dir}


# --------- Define interventions
def add_nets(cb, num_years, coverage=0.6, start_day=180):

    for year in range(num_years):
        if year % 3 == 0:
            add_ITN_age_season(cb, start=year * 365 + start_day, demographic_coverage=coverage, insecticide="bednet",
                               discard_times={'Expiration_Period_Distribution': "EXPONENTIAL_DISTRIBUTION",
                                              'Expiration_Period_Exponential': 639},
                               blocking_config=WaningEffectExponential(Decay_Time_Constant=730, Initial_Effect=0.6),
                               killing_config=WaningEffectExponential(Decay_Time_Constant=1460, Initial_Effect=0.7)
                               )

    return {'itn_coverage': coverage, 'itn_start_day': start_day}


def find_n_largest(n, demo_file_dir):

    config = configparser.ConfigParser()
    config.read('simtools.ini')
    direc = config['HPC']['input_root']

    file = os.path.join(direc, demo_file_dir, 'Demographics', 'demographics.json')
    nodes = []
    pop = []
    with open(file) as f:
        data = json.load(f)
        for node in data['Nodes']:
            nodes.append(node['NodeID'])
            pop.append(node['NodeAttributes']['InitialPopulation'])
    df = pd.DataFrame({'nodes': nodes, 'pop': pop})
    df = df.sort_values(by=['pop'], ascending=False)

    node_list = list(df['nodes'][:n])

    return node_list


def add_release(cb, demo_file_dir, number=100, num_nodes=1, release_day=180, drive_type="classic"):

    if drive_type == "classic":
        released_genome = [["X", "Y"], ["a1", "a1"]]
    elif drive_type == "integral":
        released_genome = [["X", "Y"], ["a1", "a1"], ["b1", "b1"]]

    nodelist = find_n_largest(num_nodes, demo_file_dir)

    add_mosquito_release(cb, start_day=release_day, species="gambiae", repetitions=1, number=number,
                         released_genome=released_genome, nodeIDs=nodelist)

    return {'release_number': number, 'num_nodes': num_nodes, 'release_day': release_day}


# --------- Define classic gene drive functions
# - Initial resistance frequency (defined in set_classic_genes)
# rr0 = 0  # initial frequency of resistance
# - Gene drive params (defined in add_classic_gene_drives)
# d = 0.99  # transmission rate of drive
# u = 0.5  # probability of resistance (functional)
# lne = 3 * 10 **-4  # probability of nuclease and effector loss during homing
# - Fitness params (defined in add_classic_fitness_costs)
# sd = 0  # cost of target site disruption
# sn = 0.05  # cost of nuclease expression
# se = 0.1  # cost of effector expression
# hd = 0.5  # dominance coefficient for target site disruption
# hn = 0.5  # dominance coefficient for nuclease expression
# he = 0.5  # dominance coefficient for effector expression
# - Refractoriness (defined in add_classic_fitness_costs)
# hrc = 1  # dominance coefficient for refractoriness
# rc = 1  # homozygous degree of refractorines

def set_classic_genes(cb, rr0=0):

    genes = {'gambiae': [
        {
            'Alleles': {
                'a0': 1 - rr0,  # wild type
                'a1': 0,  # complete drive construct
                'a2': rr0  # resistance allele
            },
            'Mutations': {}
        }
    ]}
    set_species_genes(cb, genes)

    return {'rr0': rr0}


def add_classic_trait_modifiers(cb, sd=0, hd=0.5, sne=0.1, hne=0.5, rc=1, hrc=1):

    traits = {'gambiae': [
        {
            'Allele_Combinations': [["a0", "a0"]],  # wt/wt
            'Trait_Modifiers': {'MORTALITY': 1}
        },
        {
            'Allele_Combinations': [["a0", "a1"]],  # wt/complete construct
            'Trait_Modifiers': {'MORTALITY': 1 / ((1 - hd * sd) * (1 - hne * sne)),
                                'INFECTED_BY_HUMAN': 1 - hrc * rc}
        },
        {
            'Allele_Combinations': [["a0", "a2"]],  # wt/resistance
            'Trait_Modifiers': {'MORTALITY': 1 / (1 - hd * sd)}
        },
        {
            'Allele_Combinations': [["a1", "a1"]],  # c/c
            'Trait_Modifiers': {'MORTALITY': 1 / ((1 - sd) * (1 - sne)),
                                'INFECTED_BY_HUMAN': 1 - rc}
        },
        {
            'Allele_Combinations': [["a1", "a2"]],  # c/r
            'Trait_Modifiers': {'MORTALITY': 1 / ((1 - sd) * (1 - hne * sne)),
                                'INFECTED_BY_HUMAN': 1 - hrc * rc}
        },
        {
            'Allele_Combinations': [["a2", "a2"]],  # r/r
            'Trait_Modifiers': {'MORTALITY': 1 / (1 - sd)}
        }
    ]}

    set_species_trait_modifiers(cb, traits)

    return {'sd': sd, 'hd': hd, 'sne': sne, 'hne': hne, 'hrc': hrc, 'rc': rc}


def add_classic_drivers(cb, d=0.99, u=0.5, lne=3e-4):

    drivers = {'gambiae': [
        {'Driver_Type': "CLASSIC",
         'Driving_Allele': "a1",
         'Alleles_Driven': [
             {
                 'Allele_To_Copy': "a1",
                 'Allele_To_Replace': "a0",
                 'Copy_To_Likelihood': {
                     'a0': (1 - d) * (1 - u),  # wild type
                     'a1': d * (1 - lne),  # complete drive construct
                     'a2': d * lne + (1 - d) * u  # resistant allele
                 }
             }
         ]
         }
    ]}

    set_species_drivers(cb, drivers)

    return {'d': d, 'u': u, 'lne': lne}


# --------- Define integral gene drive functions
# - Initial resistance frequencies (defined in set_integral_genes)
# rr10 = 0  # initial frequency of resistance (locus 1)
# rr20 = 0  # initial frequency of resistance (locus 2)
# - Gene drive params (defined in add_integral_gene_drives)
# d1 = 0.99  # transmission rate of drive
# p_nhej = 0.5  # prob of NHEJ at each locus, given no drive
# p_r_nhej = 1 / 3  # prob resistance arising from NHEJ at each locus
# p_ihdr = 1e-4  # prob of incomplete HDR at each locus, given drive
# p_r_ihdr = 1 / 3  # prob of resistance arising from incomplete HDR at each locus
# - Fitness params (defined in add_integral_fitness_costs)
# sd1 = 0  # cost of hijacking target locus 1 (drive)
# hd1 = 0.5  # dominance coefficient for hijacking at locus 1
# sd2 = 0  # cost of hijacking target locus 2 (effector)
# hd2 = 0.5  # dominance coefficient for hijacking at locus 2
# sn = 0.05  # cost of expressing nuclease at locus 1
# hn = 0.5  # dominance coefficient for expressing nuclease
# se2 = 0.1  # cost of expressing effector at locus 2
# he2 = 0.5  # dominance coefficient for expressing effector (locus 2)
# sm = 1  # cost of loss of gene function
# hm = 0.2  # dominance coefficient for loss of gene function
# - Refractoriness (defined in add_integral_fitness_costs)
# rc = 1  # homozygous degree of refractoriness
# hrc1 = 1  # dominance coefficient for refractoriness (one effector allele)

def set_integral_genes(cb, rr10=0, rr20=0):

    genes = {'gambiae': [
        {
            'Alleles': {
                'a0': 1 - rr10,
                'a1': 0,
                'a2': rr10,
                'a3': 0
            },
            'Mutations': {}
        },
        {
            'Alleles': {
                'b0': 1 - rr20,
                'b1': 0,
                'b2': rr20,
                'b3': 0
            },
            'Mutations': {}
        }
    ]}

    set_species_genes(cb, genes)

    return {'rr10': rr10, 'rr20': rr20}


def add_integral_trait_modifiers(cb, rc=1, hrc1=1, sd1=0, hd1=0.5, sd2=0, hd2=0.5,
                                 sn=0.05, hn=0.5, se2=0.1, he2=0.5, sm=1, hm=0.2):

    traits = {'gambiae': [
        {
            'Allele_Combinations': [["a0", "a0"]],
            'Trait_Modifiers': {'MORTALITY': 1}
        },
        {
            'Allele_Combinations': [["a0", "a1"]],
            'Trait_Modifiers': {'MORTALITY': 1 / ((1 - hd1 * sd1) * (1 - hn * sn))}
        },
        {
            'Allele_Combinations': [["a0", "a2"]],
            'Trait_Modifiers': {'MORTALITY': 1}
        },
        {
            'Allele_Combinations': [["a0", "a3"]],
            'Trait_Modifiers': {'MORTALITY': 1 / (1 - hm * sm)}
        },
        {
            'Allele_Combinations': [["a1", "a1"]],
            'Trait_Modifiers': {'MORTALITY': 1 / ((1 - sd1) * (1 - sn))}
        },
        {
            'Allele_Combinations': [["a1", "a2"]],
            'Trait_Modifiers': {'MORTALITY': 1 / ((1 - hd1 * sd1) * (1 - hn * sn))}
        },
        {
            'Allele_Combinations': [["a1", "a3"]],
            'Trait_Modifiers': {'MORTALITY': 1 / ((1 - hd1 * sd1) * (1 - hn * sn) * (1 - hm * sm))}
        },
        {
            'Allele_Combinations': [["a2", "a2"]],
            'Trait_Modifiers': {'MORTALITY': 1}
        },
        {
            'Allele_Combinations': [["a2", "a3"]],
            'Trait_Modifiers': {'MORTALITY': 1 / (1 - hm * sm)}
        },
        {
            'Allele_Combinations': [["a3", "a3"]],
            'Trait_Modifiers': {'ADJUST_FERTILE_EGGS': 0}
        },
        {
            'Allele_Combinations': [["b0", "b0"]],
            'Trait_Modifiers': {'MORTALITY': 1}
        },
        {
            'Allele_Combinations': [["b0", "b1"]],
            'Trait_Modifiers': {'MORTALITY': 1 / ((1 - hd2 * sd2) * (1 - he2 * se2)),
                                'INFECTED_BY_HUMAN': 1 - hrc1 * rc}
        },
        {
            'Allele_Combinations': [["b0", "b2"]],
            'Trait_Modifiers': {'MORTALITY': 1}
        },
        {
            'Allele_Combinations': [["b0", "b3"]],
            'Trait_Modifiers': {'MORTALITY': 1 / (1 - hm * sm)}
        },
        {
            'Allele_Combinations': [["b1", "b1"]],
            'Trait_Modifiers': {'MORTALITY': 1 / ((1 - sd2) * (1 - se2)),
                                'INFECTED_BY_HUMAN': 1 - rc}
        },
        {
            'Allele_Combinations': [["b1", "b2"]],
            'Trait_Modifiers': {'MORTALITY': 1 / ((1 - hd2 * sd2) * (1 - he2 * se2)),
                                'INFECTED_BY_HUMAN': 1 - hrc1 * rc}
        },
        {
            'Allele_Combinations': [["b1", "b3"]],
            'Trait_Modifiers': {
                'MORTALITY': 1 / ((1 - hd2 * sd2) * (1 - he2 * se2) * (1 - hm * sm)),
                'INFECTED_BY_HUMAN': 1 - hrc1 * rc}
        },
        {
            'Allele_Combinations': [["b2", "b2"]],
            'Trait_Modifiers': {'MORTALITY': 1}
        },
        {
            'Allele_Combinations': [["b2", "b3"]],
            'Trait_Modifiers': {'MORTALITY': 1 / (1 - hm * sm)}
        },
        {
            'Allele_Combinations': [["b3", "b3"]],
            'Trait_Modifiers': {'ADJUST_FERTILE_EGGS': 0}
        }
    ]}

    set_species_trait_modifiers(cb, traits)

    return {'rc': rc, 'hrc1': hrc1, 'sd1': sd1, 'hd1': hd1, 'sd2': sd2, 'hd2': hd2,
            'sn': sn, 'hn': hn, 'se2': se2, 'he2': he2, 'sm': sm, 'hm': hm}


def add_integral_drivers(cb, d1=0.99, p_nhej=0.5, p_ihdr=1e-4,
                         p_r_nhej=1 / 3, p_r_ihdr=1 / 3):

    drivers = {'gambiae': [{
        'Driver_Type': "INTEGRAL_AUTONOMOUS",
        'Driving_Allele': "a1",
        'Alleles_Driven': [
            {
                'Allele_To_Copy': "a1",
                'Allele_To_Replace': "a0",
                'Copy_To_Likelihood': {
                    'a0': (1 - d1) * (1 - p_nhej),
                    'a1': d1 * (1 - p_ihdr),
                    'a2': (1 - d1) * p_r_nhej * p_nhej + d1 * p_r_ihdr * p_ihdr,
                    'a3': (1 - d1) * (1 - p_r_nhej) * p_nhej + d1 * (1 - p_r_ihdr) * p_ihdr
                }
            },
            {
                'Allele_To_Copy': "b1",
                'Allele_To_Replace': "b0",
                'Copy_To_Likelihood': {
                    'b0': (1 - d1) * (1 - p_nhej),
                    'b1': d1 * (1 - p_ihdr),
                    'b2': (1 - d1) * p_r_nhej * p_nhej + d1 * p_r_ihdr * p_ihdr,
                    'b3': (1 - d1) * (1 - p_r_nhej) * p_nhej + d1 * (1 - p_r_ihdr) * p_ihdr
                }
            }
        ]
    }]
    }

    set_species_drivers(cb, drivers)

    return {'d1': d1, 'p_nhej': p_nhej, 'p_ihdr': p_ihdr,
            'p_r_nhej': p_r_nhej, 'p_r_ihdr': p_r_ihdr}
