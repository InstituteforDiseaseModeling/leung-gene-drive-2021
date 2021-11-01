import os
import numpy as np
import pandas as pd
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser

column_types = dict(Time='uint16',
                    time='uint16',
                    Run_Number='uint8',
                    Node='uint8',
                    node='uint8',
                    NodeID='uint8'
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


class AlleleFreqAnalyzer(BaseAnalyzer):

    def __init__(self, exp_name, sweep_variables, working_dir='.'):
        super(AlleleFreqAnalyzer, self).__init__(
            working_dir=working_dir,
            filenames=['output/ReportVectorGenetics_gambiae_Female_ALLELE_FREQ.csv']
        )
        self.exp_name = exp_name
        self.sweep_variables = sweep_variables
        self.num_loci = 3  # including sex chromosome
        if self.num_loci == 2:
            self.channels = ['a0', 'a1', 'a2',
                             'Total Female Population']
        elif self.num_loci == 3:
            self.channels = ['a0', 'a1', 'a2', 'a3', 'b0', 'b1', 'b2', 'b3',
                             'Total Female Population']
        self.reporter_sex = 'female'
        self.output_type = 'meanstd'  # choose: sims, meanstd
        if self.output_type == 'sims':
            self.output_fname = os.path.join(self.working_dir, "%s_spatial_avg_allele_freqs_indiv_sims" % self.exp_name)
        elif self.output_type == 'meanstd':
            self.output_fname = os.path.join(self.working_dir, "%s_spatial_avg_allele_freqs" % self.exp_name)

    def select_simulation_data(self, data, simulation):
        simdata = []

        if self.reporter_sex == 'female':
            datatemp = data['output/ReportVectorGenetics_gambiae_Female_ALLELE_FREQ.csv']
            datatemp = datatemp[datatemp['Alleles'] != 'Y']
            alleles = datatemp['Alleles'].unique()
            alleles = alleles[alleles != 'X']
            datatemp = datatemp.pivot_table('VectorPopulation', ['Time', 'NodeID'], 'Alleles').reset_index()
            datatemp = datatemp.groupby('Time').sum().reset_index()
            pop_col_name = 'Total Female Population'
            datatemp[pop_col_name] = datatemp['X'] / 2
            datatemp.drop(columns=['NodeID', 'X'], inplace=True)
            for allele in alleles:
                datatemp[allele] = datatemp[allele] / (2 * datatemp[pop_col_name])

        simdata.append(datatemp)
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
            d.to_csv(self.output_fname + fname_suffix + '.csv', index=False)

        elif self.output_type == 'meanstd':
            d_mean = d.groupby(self.sweep_variables + ['Time'])[self.channels].apply(
                np.mean).reset_index()
            d_std = d.groupby(self.sweep_variables + ['Time'])[self.channels].apply(
                np.std).reset_index()
            for channel in self.channels:
                d_mean[channel + '_std'] = d_std[channel]
            d_mean.to_csv(self.output_fname + fname_suffix + '.csv', index=False)


if __name__ == "__main__":

    SetupParser.default_block = 'HPC'
    SetupParser.init()

    out_dir = 'C:\\Users\\sleung\\OneDrive - Institute for Disease Modeling\\DTK_input_output_staging\\analyzed'

    experiments = {
        "spatialinside_vector_genetics_classic3allele_VC_and_GM_vmignegexpovnnb05_sweep_rc_d_rr0_sne":
            "cbc92044-45be-eb11-a9ec-b88303911bc1"
    }

    sweep_vars = ['rc', 'd', 'rr0', 'sne']

    for expt_name, exp_id in experiments.items():
        am = AnalyzeManager(exp_list=exp_id,
                            analyzers=[
                                AlleleFreqAnalyzer(exp_name='spatialinside_vector_genetics_classic3allele_VC_and_GM_vmignegexpovnnb05_sweep_rc_d_rr0_sne',
                                                   working_dir=out_dir,
                                                   sweep_variables=sweep_vars
                                                   )
                            ],
                            force_analyze=False)

        print(am.experiments)
        am.analyze()
