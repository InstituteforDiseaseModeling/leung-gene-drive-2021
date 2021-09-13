from dtk.utils.reports import BaseVectorGeneticsReport, BaseVectorMigrationReport, BaseVectorStatsReport
from malaria.interventions.health_seeking import add_health_seeking
from malaria.reports.MalariaReport import add_filtered_spatial_report
from malaria.study_sites.site_setup_functions import summary_report_fn
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser

from configure_interventions import *
from helper_functions import *

if __name__ == "__main__":

    ######################################
    # Set exp name
    ######################################
    # - more follow up gene drive paper runs for completion
    # exp_name = "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newsne"
    # exp_name = "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newsne"
    # exp_name = "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newse2"
    # exp_name = "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newse2"
    # exp_name = "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newsne"
    # exp_name = "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newsne"
    # exp_name = "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newse2"
    # exp_name = "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newse2"
    # exp_name = "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newsne"
    # exp_name = "spatialinside_classic3allele_GM_only_aEIR80_sweep_rc_d_rr0_sne_newsne"
    # exp_name = "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newse2"
    # exp_name = "spatialinside_integral2l4a_GM_only_aEIR80_sweep_rc_d1_rr20_se2_newse2"
    # exp_name = "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newrr0"
    # exp_name = "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newrr0"
    # exp_name = "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newrr20"
    # exp_name = "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newrr20"
    exp_name = "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newrr0"
    # exp_name = "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newrr20"

    # - follow up gene drive paper runs
    # exp_name = "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newrr0"
    # exp_name = "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newrr0"
    # exp_name = "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newrr20"
    # exp_name = "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newrr20"
    # exp_name = "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_release_number"
    # exp_name = "spatialinside_classic3allele_GM_only_aEIR30_sweep_release_number"
    # exp_name = "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_release_number"
    # exp_name = "spatialinside_integral2l4a_GM_only_aEIR30_sweep_release_number"
    # exp_name = "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_sne_newsne"
    # exp_name = "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_sne_newsne"
    # exp_name = "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_se2_newse2"
    # exp_name = "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_se2_newse2"

    # - sweeping over genetic params in diff transmission regimes, spatialinside
    # exp_name = "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_release_node_num"
    # exp_name = "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne"
    # exp_name = "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_release_day_release_node_num"
    # exp_name = "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_release_day_release_node_num"
    # exp_name = "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne"
    # exp_name = "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_release_day_release_node_num"
    # exp_name = "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2"
    # exp_name = "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2"
    # exp_name = "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2"
    # exp_name = "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2"
    # exp_name = "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2"

    # - sweeping over genetic params in diff transmission regimes, single node
    # exp_name = "single_node_vector_genetics_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_release_day_number"
    # exp_name = "single_node_vector_genetics_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_number"
    # exp_name = "single_node_vector_genetics_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_release_day_number"
    # exp_name = "single_node_vector_genetics_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_release_day_number"
    # exp_name = "single_node_vector_genetics_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_release_day_number"

    ######################################
    # Set run options
    ######################################
    run_platform = "HPC"  # choose: HPC, LOCAL
    comps_platform = "Calculon"  # choose: Calculon, Belegost
    num_cores = 15
    num_seeds = 20
    num_years = 8

    geography = "Burkina Faso"
    spatial_type = "spatial"  # choose: single_node, spatial
    input_file_dir = "BurkinaFasoSpatialInsideNodes"
    hmig_dir = "gm_def"

    exp_type = "VC_and_GM"  # choose: GM_only, VC_only, VC_and_GM, no_interventions (except health seeking)
    drive_type = "classic"  # choose: classic, integral, none
    genes_type = "classic"  # choose: classic, integral, none

    # TO ADD IN VECTOR MIGRATION INTO/OUT OF SIM REGION:
    # 1.) Set Vector_Sugar_Feeding_Frequency to VECTOR_SUGAR_FEEDING_EVERYDAY
    # in standard_cb_updates from configure_interventions
    # 2.) Uncomment add_daily_vmigin_release and add_daily_vmigout_ATSB
    # in relevant section below

    ######################################
    # Setup config file and serialized file
    ######################################
    cb = configure_VC_GM_intervention_system(geography, input_file_dir, hmig_dir,
                                             num_cores=num_cores, num_years=num_years)

    if genes_type == "classic" or genes_type == "none":
        if spatial_type == "single_node":
            # - single_node_classic3allele_aeir10_constbiting_sweep_rr0_serialization
            # serialized_exp_id = ['82c07050-4bd8-eb11-a9ec-b88303911bc1']
            # - single_node_classic3allele_aeir30_constbiting_sweep_rr0_serialization
            serialized_exp_id = ['3d9e74ad-4bd8-eb11-a9ec-b88303911bc1']
            # - single_node_classic3allele_aeir80_constbiting_sweep_rr0_serialization
            # serialized_exp_id = ['36a7b1c0-fbd5-eb11-a9ec-b88303911bc1']
        elif spatial_type == "spatial":
            # - spatialinside_classic3allele_aeir10_sweep_rr0_serialization
            # serialized_exp_id = ['752e0f20-cde8-eb11-a9ec-b88303911bc1']  # max larvhab = pow(10, 8.6)/11, 5 cores, windows/belegost
            # serialized_exp_id = ['6c7ffe6d-5904-ec11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/11, 15 cores, linux/calculon
            # serialized_exp_id = ['6086e4f3-f811-ec11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/11, oldrr0 + new rr0, 15 cores, linux/calculon
            # - spatialinside_classic3allele_aeir30_sweep_rr0_serialization
            # serialized_exp_id = ['987dc107-cde8-eb11-a9ec-b88303911bc1']  # max larvhab = pow(10, 8.6)/4.5, 5 cores, windows/belegost
            # serialized_exp_id = ['f108c2a4-76f9-eb11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/4.5, 15 cores, windows/belegost
            # serialized_exp_id = ['c39ac742-5904-ec11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/4.5, 15 cores, linux/calculon
            # serialized_exp_id = ['0deadfae-330a-ec11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/4.5, newrr0, 15 cores, linux/calculon
            # serialized_exp_id = ['16d77f12-f911-ec11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/4.5, oldrr0 + newrr0, 15 cores, linux/calculon
            # - spatialinside_classic3allele_aeir80_sweep_rr0_serialization
            # serialized_exp_id = ['3d00fcb3-bcdf-eb11-a9ec-b88303911bc1']  # max larvhab = pow(10, 8.6)/2, 5 cores, windows/belegost
            serialized_exp_id = ['ba15bf3e-f911-ec11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/2, oldrr0 + newrr0, 15 cores, linux/calculon
    elif genes_type == "integral":
        if spatial_type == "single_node":
            # - single_node_integraldrive2l4a_constbiting_sweep_rr10_rr20_serialization
            serialized_exp_id = ['7ca49772-a6cf-eb11-a9ec-b88303911bc1']
        elif spatial_type == "spatial":
            # - spatialinside_integral2l4a_aeir10_sweep_rr10_rr20_serialization
            # serialized_exp_id = ['21dae42b-e701-ec11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/11, 15 cores, linux/calculon
            # serialized_exp_id = ['2b699ff7-f80a-ec11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/11, oldrr1/20 + new rr1/20, 15 cores, linux/calculon
            # - spatialinside_integral2l4a_aeir30_sweep_rr10_rr20_serialization
            # serialized_exp_id = ['3cb08063-8efc-eb11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/4.5, 15 cores, windows/belegost
            # serialized_exp_id = ['ac84eece-fe01-ec11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/4.5, 15 cores, linux/calculon
            # serialized_exp_id = ['b9a133f5-330a-ec11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/4.5, oldrr1/20 + newrr1/20, 15 cores, linux/calculon
            # - spatialinside_integral2l4a_aeir80_sweep_rr10_rr20_serialization
            # serialized_exp_id = ['ed7673bb-e701-ec11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/2, 15 cores, linux/calculon
            serialized_exp_id = ['bf2fcc45-f90a-ec11-a9ed-b88303911bc1']  # max larvhab = pow(10, 8.6)/2, oldrr1/20 + new rr1/20, 15 cores, linux/calculon

    sim_map = []
    for exp in serialized_exp_id:
        temp_sim_map = make_simmap(exp)
        if not sim_map:
            sim_map = temp_sim_map
        else:
            sim_map = pd.concat([sim_map, temp_sim_map])

    if spatial_type == "single_node":
        serialized_file_list = ["state-14600.dtk"]
    elif spatial_type == "spatial":
        # serialized_file_list = ["state-18250-%03d.dtk" % i for i in range(0, num_cores)]
        serialized_file_list = ["state-14600-%03d.dtk" % i for i in range(0, num_cores)]

    cb.update_params({

        # - Serialization params
        'Serialization_Precision': "REDUCED",
        'Serialized_Population_Reading_Type': "READ",  # comment for no serialization
        # 'Serialized_Population_Reading_Type': "NONE",  # uncomment for no serialization
        'Serialized_Population_Writing_Type': "NONE",
        'Serialization_Mask_Node_Read': 0,
        'Enable_Random_Generator_From_Serialized_Population': 0,
        'Serialized_Population_Filenames': serialized_file_list,  # comment for no serialization

        # - Logging
        'logLevel_default': "ERROR",
        'logLevel_JsonConfigurable': "ERROR",
        'logLevel_VectorHabitat': "ERROR",
        'logLevel_Simulation': "ERROR",
        'logLevel_StandardEventCoordinator': "ERROR",
        'logLevel_LarvalHabitatMultiplier': "ERROR",
        'logLevel_SimulationEventContext': "ERROR",
        'logLevel_SusceptibilityMalaria': "ERROR"
    })

    ######################################
    # Set up sweep vars
    ######################################
    if drive_type == "classic":

        if exp_type == "VC_and_GM":

            if (exp_name == "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_release_day_release_node_num") \
                    or (exp_name == "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_release_day_release_node_num"):
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [100]
                release_node_nums = [6, 12]
                release_days = [180, 240, 300, 360, 420, 480, 545]
                rr0s = [0, 0.1, 0.2]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sds = [0]
                snes = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif (exp_name == "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne") \
                    or (exp_name == "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne"):
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr0s = [0, 0.1, 0.2]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sds = [0]
                snes = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif exp_name == "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newrr0":
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr0s = [0.001, 0.01]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sds = [0]
                snes = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif exp_name == "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_release_number":
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [1000, 10000]
                release_node_nums = [6]
                release_days = [180]
                rr0s = [0.01]
                copy_to_likelihoods = [0.95]
                sds = [0]
                snes = [0.1]
                infect_by_human_red_fracs = [0.7, 0.6]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif exp_name == "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_sne_newsne":
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr0s = [0.01]
                copy_to_likelihoods = [0.95]
                sds = [0]
                snes = [0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif (exp_name == "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newsne") \
                    or (exp_name == "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newsne") \
                    or (exp_name == "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newsne"):
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr0s = [0, 0.001, 0.01, 0.1]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sds = [0]
                snes = [0.3, 0.4, 0.5]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif (exp_name == "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newrr0") \
                    or (exp_name == "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newrr0"):
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr0s = [0.001, 0.01]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sds = [0]
                snes = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif (exp_name == "single_node_vector_genetics_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_release_day_number") \
                    or (exp_name == "single_node_vector_genetics_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_release_day_number") \
                    or (exp_name == "single_node_vector_genetics_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_release_day_number"):
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [100, 1000]
                release_node_nums = [1]
                release_days = [180, 240, 300, 360, 420, 480, 545]
                rr0s = [0, 0.1, 0.2]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sds = [0]
                snes = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]

            VC_and_GM = [
                [ModFn(DTKConfigBuilder.set_param, "Run_Number", seed),
                 ModFn(DTKConfigBuilder.set_param, "Serialized_Population_Path",
                       "{path}/output".format(path=get_outpath_for_serialized_file(
                           sim_map, comps_platform=comps_platform, genes_type=genes_type, rr0=rr0, seed=0))),
                 ModFn(add_nets, num_years=num_years, coverage=net_coverage, start_day=itn_start_day),
                 ModFn(add_release, demo_file_dir=input_file_dir, number=release_number,
                       num_nodes=release_node_num, release_day=release_day, drive_type=drive_type),
                 ModFn(set_classic_genes, rr0=rr0),
                 ModFn(add_classic_drivers, d=copy_to_likelihood),
                 ModFn(add_classic_trait_modifiers, sd=sd, sne=sne, rc=infect_by_human_red_frac),
                 ModFn(update_vmig_mul, vmig_mul=vmig_mul),
                 ModFn(update_vmig_dir, input_file_dir=input_file_dir, vmig_dir=vmig_dir)
                 ]
                for seed in range(num_seeds)
                for net_coverage in net_coverages
                for itn_start_day in itn_start_days
                for release_number in release_numbers
                for release_node_num in release_node_nums
                for release_day in release_days
                for rr0 in rr0s
                for copy_to_likelihood in copy_to_likelihoods
                for sd in sds
                for sne in snes
                for infect_by_human_red_frac in infect_by_human_red_fracs
                for vmig_mul in vmig_muls
                for vmig_dir in vmig_dirs
            ]

        elif exp_type == "GM_only":

            if (exp_name == "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_release_node_num") \
                    or (exp_name == "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_release_day_release_node_num"):
                release_numbers = [100]
                release_node_nums = [6, 12]
                release_days = [180, 240, 300, 360, 420, 480, 545]
                rr0s = [0, 0.1, 0.2]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sds = [0]
                snes = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif exp_name == "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newrr0":
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr0s = [0.001, 0.01]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sds = [0]
                snes = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif exp_name == "spatialinside_classic3allele_GM_only_aEIR30_sweep_release_number":
                release_numbers = [1000, 10000]
                release_node_nums = [6]
                release_days = [180]
                rr0s = [0.01]
                copy_to_likelihoods = [0.95]
                sds = [0]
                snes = [0.1]
                infect_by_human_red_fracs = [0.8]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif exp_name == "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_sne_newsne":
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr0s = [0.01]
                copy_to_likelihoods = [0.95]
                sds = [0]
                snes = [0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
                infect_by_human_red_fracs = [1, 0.9, 0.8]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif (exp_name == "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newsne") \
                    or (exp_name == "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newsne") \
                    or (exp_name == "spatialinside_classic3allele_GM_only_aEIR80_sweep_rc_d_rr0_sne_newsne"):
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr0s = [0, 0.001, 0.01, 0.1]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sds = [0]
                snes = [0.3, 0.4, 0.5]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif (exp_name == "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newrr0") \
                    or (exp_name == "spatialinside_classic3allele_GM_only_aEIR80_sweep_rc_d_rr0_sne_newrr0"):
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr0s = [0.001, 0.01]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sds = [0]
                snes = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif (exp_name == "single_node_vector_genetics_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_number") \
                    or (exp_name == "single_node_vector_genetics_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_release_day_number"):
                release_numbers = [100, 1000]
                release_node_nums = [1]
                release_days = [180, 240, 300, 360, 420, 480, 545]
                rr0s = [0, 0.1, 0.2]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sds = [0]
                snes = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]

            GM_only = [
                [ModFn(DTKConfigBuilder.set_param, "Run_Number", seed),
                 ModFn(DTKConfigBuilder.set_param, "Serialized_Population_Path",
                       "{path}/output".format(path=get_outpath_for_serialized_file(
                           sim_map, comps_platform=comps_platform, genes_type=genes_type, rr0=rr0, seed=0))),
                 ModFn(add_release, demo_file_dir=input_file_dir, number=release_number,
                       num_nodes=release_node_num, release_day=release_day, drive_type=drive_type),
                 ModFn(set_classic_genes, rr0=rr0),
                 ModFn(add_classic_drivers, d=copy_to_likelihood),
                 ModFn(add_classic_trait_modifiers, sd=sd, sne=sne, rc=infect_by_human_red_frac),
                 ModFn(update_vmig_mul, vmig_mul=vmig_mul),
                 ModFn(update_vmig_dir, input_file_dir=input_file_dir, vmig_dir=vmig_dir)
                 ]
                for seed in range(num_seeds)
                for release_number in release_numbers
                for release_node_num in release_node_nums
                for release_day in release_days
                for rr0 in rr0s
                for copy_to_likelihood in copy_to_likelihoods
                for sd in sds
                for sne in snes
                for infect_by_human_red_frac in infect_by_human_red_fracs
                for vmig_mul in vmig_muls
                for vmig_dir in vmig_dirs
            ]

    elif drive_type == "integral":

        if exp_type == "VC_and_GM":

            if (exp_name == "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2") \
                    or (exp_name == "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2") \
                    or (exp_name == "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2"):
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr10s = [0]
                rr20s = [0, 0.1, 0.2]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sns = [0.05]
                se2s = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif exp_name == "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_release_number":
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [1000, 10000]
                release_node_nums = [6]
                release_days = [180]
                rr10s = [0]
                rr20s = [0.01]
                copy_to_likelihoods = [0.95]
                sns = [0.05]
                se2s = [0.1]
                infect_by_human_red_fracs = [0.7, 0.6]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif exp_name == "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newrr20":
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr10s = [0]
                rr20s = [0.001, 0.01]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sns = [0.05]
                se2s = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif exp_name == "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_se2_newse2":
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr10s = [0]
                rr20s = [0.01]
                copy_to_likelihoods = [0.95]
                sns = [0.05]
                se2s = [0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif (exp_name == "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newse2") \
                    or (exp_name == "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newse2") \
                    or (exp_name == "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newse2"):
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr10s = [0]
                rr20s = [0, 0.001, 0.01, 0.1]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sns = [0.05]
                se2s = [0.3, 0.4, 0.5]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif (exp_name == "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newrr20") \
                    or (exp_name == "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newrr20"):
                net_coverages = [0.7]
                itn_start_days = [180]
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr10s = [0]
                rr20s = [0.001, 0.01]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sns = [0.05]
                se2s = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            VC_and_GM = [
                [ModFn(DTKConfigBuilder.set_param, "Run_Number", seed),
                 ModFn(DTKConfigBuilder.set_param, "Serialized_Population_Path",
                       "{path}/output".format(path=get_outpath_for_serialized_file(
                           sim_map, comps_platform=comps_platform, genes_type=genes_type, rr10=rr10, rr20=rr20, seed=0))),
                 ModFn(add_nets, num_years=num_years, coverage=net_coverage, start_day=itn_start_day),
                 ModFn(add_release, demo_file_dir=input_file_dir, number=release_number,
                       num_nodes=release_node_num, release_day=release_day, drive_type=drive_type),
                 ModFn(set_integral_genes, rr10=rr10, rr20=rr20),
                 ModFn(add_integral_drivers, d1=copy_to_likelihood),
                 ModFn(add_integral_trait_modifiers, rc=infect_by_human_red_frac, sn=sn, se2=se2),
                 ModFn(update_vmig_mul, vmig_mul=vmig_mul),
                 ModFn(update_vmig_dir, input_file_dir=input_file_dir, vmig_dir=vmig_dir)
                 ]
                for seed in range(num_seeds)
                for net_coverage in net_coverages
                for itn_start_day in itn_start_days
                for release_number in release_numbers
                for release_node_num in release_node_nums
                for release_day in release_days
                for rr10 in rr10s
                for rr20 in rr20s
                for copy_to_likelihood in copy_to_likelihoods
                for infect_by_human_red_frac in infect_by_human_red_fracs
                for sn in sns
                for se2 in se2s
                for vmig_mul in vmig_muls
                for vmig_dir in vmig_dirs
            ]

        elif exp_type == "GM_only":

            if (exp_name == "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2") \
                    or (exp_name == "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2"):
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr10s = [0]
                rr20s = [0, 0.1, 0.2]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sns = [0.05]
                se2s = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif exp_name == "spatialinside_integral2l4a_GM_only_aEIR30_sweep_release_number":
                release_numbers = [1000, 10000]
                release_node_nums = [6]
                release_days = [180]
                rr10s = [0]
                rr20s = [0.01]
                copy_to_likelihoods = [0.95]
                sns = [0.05]
                se2s = [0.1]
                infect_by_human_red_fracs = [0.8]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif exp_name == "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newrr20":
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr10s = [0]
                rr20s = [0.001, 0.01]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sns = [0.05]
                se2s = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif exp_name == "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_se2_newse2":
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr10s = [0]
                rr20s = [0.01]
                copy_to_likelihoods = [0.95]
                sns = [0.05]
                se2s = [0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
                infect_by_human_red_fracs = [1, 0.9, 0.8]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif (exp_name == "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newse2") \
                    or (exp_name == "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newse2") \
                    or (exp_name == "spatialinside_integral2l4a_GM_only_aEIR80_sweep_rc_d1_rr20_se2_newse2"):
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr10s = [0]
                rr20s = [0, 0.001, 0.01, 0.1]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sns = [0.05]
                se2s = [0.3, 0.4, 0.5]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            elif (exp_name == "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newrr20") \
                    or (exp_name == "spatialinside_integral2l4a_GM_only_aEIR80_sweep_rc_d1_rr20_se2_newrr20"):
                release_numbers = [100]
                release_node_nums = [6]
                release_days = [180]
                rr10s = [0]
                rr20s = [0.001, 0.01]
                copy_to_likelihoods = [1, 0.95, 0.9]
                sns = [0.05]
                se2s = [0, 0.05, 0.1, 0.15, 0.2]
                infect_by_human_red_fracs = [1, 0.9, 0.8, 0.7, 0.6, 0.5]
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b05"]

            GM_only = [
                [ModFn(DTKConfigBuilder.set_param, "Run_Number", seed),
                 ModFn(DTKConfigBuilder.set_param, "Serialized_Population_Path",
                       "{path}/output".format(path=get_outpath_for_serialized_file(
                           sim_map, comps_platform=comps_platform, genes_type=genes_type, rr10=rr10, rr20=rr20, seed=0))),
                 ModFn(add_release, demo_file_dir=input_file_dir, number=release_number,
                       num_nodes=release_node_num, release_day=release_day, drive_type=drive_type),
                 ModFn(set_integral_genes, rr10=rr10, rr20=rr20),
                 ModFn(add_integral_drivers, d1=copy_to_likelihood),
                 ModFn(add_integral_trait_modifiers, rc=infect_by_human_red_frac, sn=sn, se2=se2),
                 ModFn(update_vmig_mul, vmig_mul=vmig_mul),
                 ModFn(update_vmig_dir, input_file_dir=input_file_dir, vmig_dir=vmig_dir)
                 ]
                for seed in range(num_seeds)
                for release_number in release_numbers
                for release_node_num in release_node_nums
                for release_day in release_days
                for rr10 in rr10s
                for rr20 in rr20s
                for copy_to_likelihood in copy_to_likelihoods
                for infect_by_human_red_frac in infect_by_human_red_fracs
                for sn in sns
                for se2 in se2s
                for vmig_mul in vmig_muls
                for vmig_dir in vmig_dirs
            ]

    elif drive_type == "none":

        if exp_type == "VC_only":

            if genes_type == "classic":

                # - spatialinside_vector_genetics_classic3allele_VC_only_sweep_rr0_vmigdir"
                net_coverages = [0.7]
                itn_start_days = [180]
                rr0s = [0, 0.05, 0.1, 0.25, 0.4]
                vmig_muls = [1]
                vmig_dirs = ["gm_def", "thomas2013negexpovnn_b05"]

                VC_only = [
                    [ModFn(DTKConfigBuilder.set_param, "Run_Number", seed),
                     ModFn(DTKConfigBuilder.set_param, "Serialized_Population_Path",
                           "{path}/output".format(path=get_outpath_for_serialized_file(
                               sim_map, comps_platform=comps_platform, genes_type=genes_type, rr0=rr0, seed=0))),
                     ModFn(add_nets, num_years=num_years, coverage=net_coverage, start_day=itn_start_day),
                     ModFn(set_classic_genes, rr0=rr0),
                     ModFn(update_vmig_mul, vmig_mul=vmig_mul),
                     ModFn(update_vmig_dir, input_file_dir=input_file_dir, vmig_dir=vmig_dir)
                     ]
                    for seed in range(num_seeds)
                    for itn_start_day in itn_start_days
                    for net_coverage in net_coverages
                    for rr0 in rr0s
                    for vmig_mul in vmig_muls
                    for vmig_dir in vmig_dirs
                ]

            elif genes_type == "integral":

                # - spatialinside_vector_genetics_integral2l4a_VC_only_sweep_rr10_rr20_vmigdir"
                net_coverages = [0.7]
                itn_start_days = [180]
                rr10s = [0, 0.05, 0.1, 0.25, 0.4]
                rr20s = [0, 0.05, 0.1, 0.25, 0.4]
                vmig_muls = [1]
                vmig_dirs = ["gm_def", "thomas2013negexpovnn_b05"]

                VC_only = [
                    [ModFn(DTKConfigBuilder.set_param, "Run_Number", seed),
                     ModFn(DTKConfigBuilder.set_param, "Serialized_Population_Path",
                           "{path}/output".format(path=get_outpath_for_serialized_file(
                               sim_map, comps_platform=comps_platform, genes_type=genes_type, rr10=rr10, rr20=rr20, seed=0))),
                     ModFn(add_nets, num_years=num_years, coverage=net_coverage, start_day=itn_start_day),
                     ModFn(set_integral_genes, rr10=rr10, rr20=rr20),
                     ModFn(update_vmig_mul, vmig_mul=vmig_mul),
                     ModFn(update_vmig_dir, input_file_dir=input_file_dir, vmig_dir=vmig_dir)
                     ]
                    for seed in range(num_seeds)
                    for itn_start_day in itn_start_days
                    for net_coverage in net_coverages
                    for rr10 in rr10s
                    for rr20 in rr20s
                    for vmig_mul in vmig_muls
                    for vmig_dir in vmig_dirs
                ]

        elif exp_type == "no_interventions":

            if genes_type == "none":

                # - spatialinsidenodeID71LH_no_srlztn_vmigthomas2013negexpovnn_sweep_vmig_dir
                vmig_muls = [1]
                vmig_dirs = ["thomas2013negexpovnn_b03", "thomas2013negexpovnn_b04",
                             "thomas2013negexpovnn_b05", "thomas2013negexpovnn_b06",
                             "thomas2013negexpovnn_b07", "thomas2013negexpovnn_b1"]

                no_interventions = [
                    [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
                     # ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path',
                     #       "{path}/output".format(path=get_outpath_for_serialized_file(
                     #           sim_map, comps_platform=comps_platform, genes_type=genes_type, rr0=0, seed=0))),
                     ModFn(update_vmig_mul, vmig_mul=vmig_mul),
                     ModFn(update_vmig_dir, input_file_dir=input_file_dir, vmig_dir=vmig_dir),
                     ]
                    for seed in range(num_seeds)
                    for vmig_mul in vmig_muls
                    for vmig_dir in vmig_dirs
                ]

            elif genes_type == "classic":

                no_interventions = [
                    [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
                     ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path',
                           "{path}/output".format(path=get_outpath_for_serialized_file(
                               sim_map, comps_platform=comps_platform, genes_type=genes_type, rr0=rr0, seed=0))),
                     ModFn(set_classic_genes, rr0=rr0),
                     ]
                    for seed in range(num_seeds)
                    for rr0 in rr0s
                ]

            elif genes_type == "integral":

                no_interventions = [
                    [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
                     ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path',
                           "{path}/output".format(path=get_outpath_for_serialized_file(
                               sim_map, comps_platform=comps_platform, genes_type=genes_type, rr10=rr10, rr20=rr20, seed=0))),
                     ModFn(set_integral_genes, rr10=rr10, rr20=rr20),
                     ]
                    for seed in range(num_seeds)
                    for rr10 in rr10s
                    for rr20 in rr20s
                ]

    ######################################
    # Define model builder
    ######################################
    if exp_type == "no_interventions":
        builder = ModBuilder.from_list(no_interventions)
    elif exp_type == "GM_only":
        builder = ModBuilder.from_list(GM_only)
    elif exp_type == "VC_only":
        builder = ModBuilder.from_list(VC_only)
    elif exp_type == "VC_and_GM":
        builder = ModBuilder.from_list(VC_and_GM)

    ######################################
    # Set up health seeking interventions
    ######################################
    add_health_seeking(cb,
                       targets=[{'trigger': "NewClinicalCase",
                                 'coverage': 0.5,
                                 'agemin': 0,
                                 'agemax': 100,
                                 'seek': 1,
                                 'rate': 0.3},
                                {'trigger': "NewSevereCase",
                                 'coverage': 0.8,
                                 'agemin': 0,
                                 'agemax': 100,
                                 'seek': 1,
                                 'rate': 0.5}
                                ],
                       drug=["Artemether", "Lumefantrine"],
                       start_day=0,
                       broadcast_event_name="Received_Treatment")

    ######################################
    # Set up reports
    ######################################
    # ---- Summary report
    # summary_report_fn(start=0, interval=30.0,
    #                   description='Monthly_Report',
    #                   age_bins=[5.0, 15.0, 50.0, 100.0],
    #                   parasitemia_bins=[0.0, 100.0, 4000000.0],
    #                   infection_bins=np.arange(0.0, 101.0, 5.0).tolist(),
    #                   ipfilter='')(cb)

    # ---- All genomes into 1 file
    # cb.add_reports(BaseVectorGeneticsReport(type="ReportVectorGenetics",
    #                                         species="gambiae",
    #                                         gender="VECTOR_FEMALE",
    #                                         include_vector_state_columns=0,
    #                                         stratify_by="GENOME",
    #                                         combine_similar_genomes=1
    #                                         ))

    # ---- All allele frequencies into 1 file
    cb.add_reports(BaseVectorGeneticsReport(type='ReportVectorGenetics',
                                            species='gambiae',
                                            gender='VECTOR_FEMALE',
                                            # gender='VECTOR_BOTH_GENDERS',
                                            include_vector_state_columns=0,
                                            stratify_by='ALLELE_FREQ'
                                            ))

    # ---- Vector migration
    # cb.add_reports(BaseVectorMigrationReport(type='ReportVectorMigration',
    #                                          start_day=224, end_day=253
    #                                          ))

    # ---- Spatial report(s) OR can update config file params:
    # add_filtered_spatial_report(cb, start=0, end=365 * num_years,
    #                             channels=["Population", "Prevalence", "New_Clinical_Cases",
    #                                       "Daily_EIR", "Adult_Vectors"])

    ######################################
    # Run
    ######################################
    run_sim_args = {'config_builder': cb,
                    'exp_name': exp_name,
                    'exp_builder': builder}

    SetupParser.default_block = run_platform

    SetupParser.init(run_platform)
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
