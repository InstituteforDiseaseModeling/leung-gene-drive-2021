import os
import pandas as pd
import numpy as np
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser

column_types = dict(Time='uint16',
                    time='uint16',
                    Run_Number='uint8'
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


class InsetAnalyzer(BaseAnalyzer):

    def __init__(self, exp_name, report_names=["InsetChart"], channels=None, sweep_variables=None, working_dir="."):
        super(InsetAnalyzer, self).__init__(working_dir=working_dir, filenames=["output/InsetChart.json"])
        self.sweep_variables = sweep_variables
        self.channels = channels
        self.reports = report_names
        self.exp_name = exp_name
        self.output_type = 'meanstd'  # choose: sims, meanstd, both
        self.output_fname_sims = os.path.join(self.working_dir, "%s_inset_data_indiv_sims" % self.exp_name)
        self.output_fname_meanstd = os.path.join(self.working_dir, "%s_inset_data" % self.exp_name)

    def select_simulation_data(self, data, simulation):
        simdata = []

        datatemp = data["output/InsetChart.json"]

        adult_vectors = pd.Series(datatemp['Channels']['Adult Vectors']['Data'])
        infectious_vectors = pd.Series(datatemp['Channels']['Infectious Vectors']['Data'])
        hrp2_prevalence = pd.Series(datatemp['Channels']['PfHRP2 Prevalence']['Data'])
        true_prevalence = pd.Series(datatemp['Channels']['True Prevalence']['Data'])

        df = pd.DataFrame(columns=[
            'Adult Vectors', 'Infectious Vectors',
            'PfHRP2 Prevalence', 'True Prevalence',
        ])
        df['Adult Vectors'] = adult_vectors
        df['Infectious Vectors'] = infectious_vectors
        df['PfHRP2 Prevalence'] = hrp2_prevalence
        df['True Prevalence'] = true_prevalence

        df['Time'] = df.index

        simdata.append(df)
        simdata = pd.concat(simdata)

        for sweep_var in self.sweep_variables:
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

        if self.output_type == 'sims':
            d.to_csv(self.output_fname_sims + fname_suffix + '.csv', index=False)

        elif self.output_type == 'meanstd':
            d_mean = d.groupby(self.sweep_variables + ['Time'])[self.channels].apply(
                np.mean).reset_index()
            d_std = d.groupby(self.sweep_variables + ['Time'])[self.channels].apply(
                np.std).reset_index()
            for channel in self.channels:
                d_mean[channel + '_std'] = d_std[channel]
            d_mean.to_csv(self.output_fname_meanstd + fname_suffix + '.csv', index=False)

        elif self.output_type == 'both':
            d.to_csv(self.output_fname_sims + fname_suffix + '.csv', index=False)
            d_mean = d.groupby(self.sweep_variables + ['Time'])[self.channels].apply(np.mean).reset_index()
            d_std = d.groupby(self.sweep_variables + ['Time'])[self.channels].apply(np.std).reset_index()
            for channel in self.channels:
                d_mean[channel + '_std'] = d_std[channel]
            d_mean.to_csv(self.output_fname_meanstd + fname_suffix + '.csv', index=False)


if __name__ == "__main__":

    SetupParser.default_block = 'HPC'
    SetupParser.init()

    experiments = {
        "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne":
            "69852bb5-0b05-ec11-a9ed-b88303911bc1"
    }

    sweep_vars = ['rc', 'd', 'rr0', 'sne']

    for exp_name, exp_id in experiments.items():
        am = AnalyzeManager(exp_list=exp_id,
                            analyzers=[
                                InsetAnalyzer(exp_name=exp_name,
                                              sweep_variables=sweep_vars,
                                              channels=['Adult Vectors', 'Infectious Vectors',
                                                        'True Prevalence', 'PfHRP2 Prevalence']
                                              )
                                       ],
                            force_analyze=False)

        print(am.experiments)
        am.analyze()

