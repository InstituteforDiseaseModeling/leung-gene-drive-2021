import os

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_species, set_species_param


def standard_cb_updates(cb, geography, input_file_dir, hmig_dir):

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

    # - Climate
    cb.update_params({
        'Climate_Model': "CLIMATE_CONSTANT",
    })

    # - Necessary specs
    cb.update_params({
        'Inset_Chart_Reporting_Include_30Day_Avg_Infection_Duration': 1,
        'Enable_Malaria_CoTransmission': 0
    })


def update_vector_params(cb):

    cb.update_params({'Vector_Species_Names': ["gambiae"]})
    set_species(cb, ["gambiae"])
    set_species_param(cb, "gambiae", "Anthropophily", 0.65)

    Insecticides = [
        {
            "Name": "bednet",
            "Resistances": []
        }
    ]

    cb.update_params({
        'Egg_Hatch_Density_Dependence': "NO_DENSITY_DEPENDENCE",
        'Temperature_Dependent_Feeding_Cycle': "NO_TEMPERATURE_DEPENDENCE",
        'Enable_Drought_Egg_Hatch_Delay': 0,
        'Enable_Egg_Mortality': 0,
        'Enable_Temperature_Dependent_Egg_Hatching': 0,
        'Insecticides': Insecticides,
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
                      }})
    set_species_param(cb, "gambiae", "Adult_Life_Expectancy", 20)
    set_species_param(cb, "gambiae", "Male_Life_Expectancy", 10)
    set_species_param(cb, "gambiae", "Indoor_Feeding_Fraction", 0.9)
    set_species_param(cb, "gambiae", "Vector_Sugar_Feeding_Frequency",
                      "VECTOR_SUGAR_FEEDING_NONE")


def configure_VC_GM_intervention_system(geography, input_file_dir, hmig_dir,
                                        num_cores=1, num_years=100):

    cb = DTKConfigBuilder.from_defaults("MALARIA_SIM",
                                        Num_Cores=num_cores,
                                        Simulation_Duration=int(365 * num_years)
                                        )
    standard_cb_updates(cb, geography, input_file_dir, hmig_dir)
    update_vector_params(cb)

    return cb
