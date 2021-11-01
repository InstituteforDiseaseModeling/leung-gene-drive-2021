import numpy as np
import os

from dtk.generic.serialization import add_SerializationTimesteps
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.utils.reports import BaseVectorGeneticsReport, BaseVectorMigrationReport, BaseVectorStatsReport
from dtk.vector.species import set_species, set_species_param
from malaria.interventions.health_seeking import add_health_seeking
from malaria.study_sites.site_setup_functions import summary_report_fn
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser

from helper_functions import set_classic_genes, set_integral_genes, update_vmig_mul

if __name__ == "__main__":

    ######################################
    # Set exp name
    ######################################
    # - spatial, varying eir
    # exp_name = "spatialinside_integral2l4a_aeir30_sweep_rr10_rr20_serialization"  # max larvhab=pow(10, 8.6)/4.5
    # exp_name = "spatialinside_integral2l4a_aeir10_sweep_rr10_rr20_serialization"  # max larvhab=pow(10, 8.6)/11
    # exp_name = "spatialinside_integral2l4a_aeir80_sweep_rr10_rr20_serialization"  # max larvhab=pow(10, 8.6)/2
    # exp_name = "spatialinside_classic3allele_aeir30_sweep_rr0_serialization"  # max larvhab=pow(10, 8.6)/5, /4 2nd time, /4.5 3rd time
    # exp_name = "spatialinside_classic3allele_aeir10_sweep_rr0_serialization"  # max larvhab=pow(10, 8.6)/15, /8 2nd time, /11 3rd time
    exp_name = "spatialinside_classic3allele_aeir80_sweep_rr0_serialization"  # max larvhab=pow(10, 8.6)/2

    # - single node, varying eir
    # exp_name = "single_node_classic3allele_aeir30_constbiting_sweep_rr0_serialization" # max larvhab=pow(10, 8.6)/5
    # exp_name = "single_node_classic3allele_aeir10_constbiting_sweep_rr0_serialization" # max larvhab=pow(10, 8.6)/16
    # exp_name = "single_node_classic3allele_aeir80_constbiting_sweep_rr0_serialization" # max larvhab=pow(10, 8.6)/2

    ######################################
    # Set run options
    ######################################
    run_platform = "HPC"  # choose: HPC, LOCAL

    # (should be less than or equal to number of nodes requested by
    # summary report and should ensure each core gets a job)
    num_cores = 15
    num_seeds = 1
    num_years = 40

    geography = "Burkina Faso"
    input_file_dir = "BurkinaFasoSpatialInsideNodes"  # (check simtools.ini for full path; this is the highest dir)
    hmig_dir = "gm_def"
    vmig_dir = "thomas2013negexpovnn_b05"

    genes_type = "classic"  # choose: classic, integral

    ######################################
    # Setup config file and serialized file
    ######################################
    cb = DTKConfigBuilder.from_defaults('MALARIA_SIM',
                                        Num_Cores=num_cores,
                                        Simulation_Duration=int(365 * num_years))

    # - Serialization params
    add_SerializationTimesteps(cb, [num_years * 365], end_at_final=True)
    cb.update_params({
        'Serialization_Precision': "REDUCED",
        'Serialized_Population_Reading_Type': "NONE",
        'Serialization_Mask_Node_Write': 0,
        'Serialized_Population_Writing_Type': "TIMESTEP",
    })

    # - Logging
    cb.update_params({
        'logLevel_JsonConfigurable': "ERROR",
        'logLevel_VectorHabitat': "ERROR",
        'logLevel_Simulation': "ERROR",
        'logLevel_StandardEventCoordinator': "ERROR",
        'logLevel_LarvalHabitatMultiplier': "ERROR",
        'logLevel_SimulationEventContext': "ERROR",
        'logLevel_SusceptibilityMalaria': "ERROR"
    })

    # - Demographics and geography
    cb.update_params({
        'Demographics_Filenames': [os.path.join(input_file_dir, "Demographics", "demographics.json")],
        'Geography': geography,
        'Age_Initialization_Distribution_Type': "DISTRIBUTION_COMPLEX",
        'Birth_Rate_Dependence': "FIXED_BIRTH_RATE",
        'Enable_Nondisease_Mortality': 1,
        'Enable_Demographics_Risk': 1,
        'New_Diagnostic_Sensitivity': 0.025,  # 40/uL
        'Default_Geography_Initial_Node_Population': 1000,
        'Default_Geography_Torus_Size': 10,
        'Disable_IP_Whitelist': 1,
        'Disable_NP_Whitelist': 1
    })
    cb.set_param("Enable_Demographics_Builtin", 0)
    cb.set_param("Valid_Intervention_States", [])

    # - Vector migration
    if vmig_dir is not None:
        cb.update_params({
            'Vector_Migration_Filename_Local': os.path.join(input_file_dir, "VectorMigration", vmig_dir,
                                                            "vector_local_migration.bin"),
            'Enable_Vector_Migration': 1,
            'Enable_Vector_Migration_Local': 1,
            'Enable_Vector_Migration_Human': 0,
            'Enable_Vector_Migration_Wind': 0,
            'Vector_Migration_Food_Modifier': 0,
            'Vector_Migration_Habitat_Modifier': 0,
            'Vector_Migration_Modifier_Equation': "LINEAR",
            'Vector_Migration_Stay_Put_Modifier': 0
        })

    # - Human migration
    if hmig_dir is not None:
        cb.update_params({
            'Roundtrip_Waypoints': 0,
            'Local_Migration_Filename': os.path.join(input_file_dir, "HumanMigration", hmig_dir, "human_local_migration.bin"),
            'Enable_Local_Migration': 1,
            'Migration_Model': "FIXED_RATE_MIGRATION",
            'Migration_Pattern': "SINGLE_ROUND_TRIPS",
            'Local_Migration_Roundtrip_Duration': 2,  # mean of exponential days-at-destination distribution
            'Local_Migration_Roundtrip_Probability': 0.95,  # fraction that return
            'x_Local_Migration': 10
        })

    # - Climate settings
    cb.update_params(({
        'Climate_Model': "CLIMATE_CONSTANT"
    }))

    # - Necessary specs
    cb.update_params({
        'Inset_Chart_Reporting_Include_30Day_Avg_Infection_Duration': 1,
        'Enable_Malaria_CoTransmission': 0
    })

    # - Vector parameters
    cb.update_params({
        'Vector_Species_Names': ["gambiae"],
        'Vector_Sampling_Type': "TRACK_ALL_VECTORS"
    })
    set_species(cb, ["gambiae"])
    set_species_param(cb, "gambiae", "Anthropophily", 0.65)

    cb.update_params({
        'Egg_Hatch_Density_Dependence': "NO_DENSITY_DEPENDENCE",
        'Temperature_Dependent_Feeding_Cycle': "NO_TEMPERATURE_DEPENDENCE",
        'Enable_Drought_Egg_Hatch_Delay': 0,
        'Enable_Egg_Mortality': 0,
        'Enable_Temperature_Dependent_Egg_Hatching': 0,
    })

    set_species_param(cb, "gambiae", "Larval_Habitat_Types",
                      {'LINEAR_SPLINE': {
                          'Capacity_Distribution_Over_Time': {
                              'Times': [0.0, 30.417, 60.833, 91.25, 121.667, 152.083,
                                        182.5, 212.917, 243.333, 273.75, 304.167, 334.583],
                              'Values': [3, 0.8, 1.25, 0.1, 2.7, 10, 6, 35, 2.8, 1.5, 1.6, 2.1]
                          },
                          'Capacity_Distribution_Number_Of_Years': 1,
                          'Max_Larval_Capacity': pow(10, 8.6)/2
                      }}
                      )
    set_species_param(cb, "gambiae", "Adult_Life_Expectancy", 20)
    set_species_param(cb, "gambiae", "Male_Life_Expectancy", 10)
    set_species_param(cb, "gambiae", "Indoor_Feeding_Fraction", 0.9)
    set_species_param(cb, "gambiae", "Vector_Sugar_Feeding_Frequency",
                      "VECTOR_SUGAR_FEEDING_NONE")

    ########################## VECTOR GENETICS ####################################################
    if genes_type == "classic":
        builder = ModBuilder.from_list(
            [
                [ModFn(DTKConfigBuilder.set_param, "Run_Number", seed),
                 ModFn(set_classic_genes, rr0=rr0),
                 ]
                for seed in range(num_seeds)
                for rr0 in [0, 0.001, 0.01, 0.1]
            ]
        )

    elif genes_type == "integral":
        builder = ModBuilder.from_list(
            [
                [ModFn(DTKConfigBuilder.set_param, "Run_Number", seed),
                 ModFn(set_integral_genes, rr10=rr10, rr20=rr20)
                 ]
                for seed in range(num_seeds)
                for rr10 in [0, 0.001, 0.01, 0.1]
                for rr20 in [0, 0.001, 0.01, 0.1]
            ]
        )

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
                       start_day=(num_years - 15) * 365,
                       broadcast_event_name="Received_Treatment")

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
