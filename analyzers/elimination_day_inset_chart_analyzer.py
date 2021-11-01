import os
import pandas as pd
import numpy as np
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser

column_types = dict(Run_Number='uint8',
                    True_Prevalence_elim='bool'
                    )


def optimize_dataframe(simdata):
    """
    Set the data to the optimal size for their types. For examples, we know this experiment does not have
    more than 255 nodes so we should assign the NodeId columns to an uint8 since it takes less memory than the
    default uint64
    :param simdata:
    :return:
    """
    for column, item_type in column_types.items():
        if column in simdata.columns:
            simdata[column] = simdata[column].astype(item_type)
    return simdata


class EliminationDayInsetAnalyzer(BaseAnalyzer):

    def __init__(self, exp_name, report_names=["InsetChart"], sweep_variables=None, working_dir="."):
        super(EliminationDayInsetAnalyzer, self).__init__(working_dir=working_dir, filenames=["output/InsetChart.json"])
        self.sweep_variables = sweep_variables
        self.reports = report_names
        self.exp_name = exp_name
        self.output_fname_sims = os.path.join(
            self.working_dir, "%s_inset_data_elim_day_number_indiv_sims"
                              % self.exp_name)

    def select_simulation_data(self, data, simulation):
        simdata = []

        datatemp = data["output/InsetChart.json"]

        true_prevalence = pd.Series(datatemp['Channels']['True Prevalence']['Data'])

        if (true_prevalence == 0).any():
            elim_day = true_prevalence[true_prevalence == 0].index[0]
            df = pd.DataFrame(
                [[
                    (true_prevalence[elim_day:] == 0).all(), elim_day
                ]],
                columns=['True_Prevalence_elim', 'True_Prevalence_elim_day'])

        else:
            df = pd.DataFrame(
                [[
                    False, np.nan
                ]],
                columns=['True_Prevalence_elim', 'True_Prevalence_elim_day'])

        simdata.append(df)
        simdata = pd.concat(simdata)

        for sweep_var in self.sweep_variables + ['Run_Number']:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
            else:
                simdata[sweep_var] = 0

        simdata.reset_index(drop=True)
        return optimize_dataframe(simdata)

    def finalize(self, d: pd.DataFrame, partition_vars: list, partition_vars_vals: list):

        fname_suffix = ''
        for ipv, partition_var in enumerate(partition_vars):
            fname_suffix = fname_suffix + '_' + partition_var + str(partition_vars_vals[ipv])

        d.to_csv(self.output_fname_sims + fname_suffix + '.csv', index=False)


if __name__ == "__main__":
    SetupParser.default_block = 'HPC'
    SetupParser.init()

    projectdir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\DTK_input_output_staging\\analyzed'
    out_dir = os.path.join(projectdir)

    experiments = {
        "inputEIRnode2_2locus_2nodes_nomigrate_sweepitncov_sweepxtlh":
            "e2c91b14-b114-eb11-a2c7-c4346bcb1553"
    }

    sweep_vars = ['Run_Number', 'ITN_coverage', 'x_Temporary_Larval_Habitat']

    for exp_name, exp_id in experiments.items():
        am = AnalyzeManager(exp_list=exp_id,
                            analyzers=[
                                EliminationDayInsetAnalyzer(working_dir=out_dir,
                                                            exp_name=exp_name,
                                                            sweep_variables=sweep_vars)
                            ],
                            force_analyze=False)

        print(am.experiments)
        am.analyze()
