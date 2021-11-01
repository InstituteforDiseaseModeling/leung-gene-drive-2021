import os
import sys

sys.path.append(os.path.dirname(__file__))

from partitioning_analyzer import PartitioningDataAnalyzeManager

from allele_freq_spatial_avg_analyzer import AlleleFreqAnalyzer
from elimination_day_inset_chart_analyzer import EliminationDayInsetAnalyzer
from elimination_set_day_inset_chart_analyzer import EliminationSetDayInsetAnalyzer
from inset_chart_analyzer import InsetAnalyzer


if __name__ == "__main__":
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne":
    #         "a661f6b1-5525-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newsne":
    #         "3bef088d-ab25-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newrr0":
    #         "8f29c9da-ab25-ec11-9ecd-9440c9bee941"
    # }

    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne":
    #         "940c6484-3626-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newsne":
    #         "1b495c58-ac25-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newrr0":
    #         "c004e0ed-ac25-ec11-9ecd-9440c9bee941"
    # }

    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2":
    #         "87952b20-b221-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newse2":
    #         "05e863ae-7f12-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newrr20":
    #         "4578e934-a312-ec11-a9ed-b88303911bc1"
    # }

    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2":
    #         "1a3f69cf-b025-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newse2":
    #         "f5ab2c5b-8012-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newrr20":
    #         "81a418c1-a212-ec11-a9ed-b88303911bc1"
    # }

    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne":
    #         "f5a059df-ae21-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newsne":
    #         "e265629a-a725-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newrr0":
    #         "253b91f5-2826-ec11-9ecd-9440c9bee941"
    # }

    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne":
    #         "4abc9dfa-2a26-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newsne":
    #         "e6a98024-a925-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newrr0":
    #         "755c1724-3626-ec11-9ecd-9440c9bee941"
    # }

    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2":
    #         "7fdeaca5-eefe-eb11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newse2":
    #         "9250cec6-af21-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newrr20":
    #         "b2bdccfc-b80a-ec11-a9ed-b88303911bc1"
    # }

    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2":
    #         "1e2dd6d7-2926-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newse2":
    #         "79f3816b-b021-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newrr20":
    #         "32d788f3-b70a-ec11-a9ed-b88303911bc1"
    # }

    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne":
    #         "0ee573f0-5625-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newsne":
    #         "65b564cf-b125-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newrr0":
    #         "bad4ac1d-b225-ec11-9ecd-9440c9bee941"
    # }

    experiments = {
        "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2":
            "2810b14b-b125-ec11-9ecd-9440c9bee941"
    }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newrr20":
    #         "e91d8c88-a312-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newse2":
    #         "b4b8bbb4-8012-ec11-a9ed-b88303911bc1"
    # }

    # sweep_vars = ['rc', 'd', 'rr0', 'sne']
    sweep_vars = ['rc', 'd1', 'rr20', 'se2']

    for expt_name, exp_id in experiments.items():
        am = PartitioningDataAnalyzeManager(exp_list=exp_id,
                                            partitionable_columns=sweep_vars,
                                            analyzers=
                                            [
                                                AlleleFreqAnalyzer(
                                                    exp_name=expt_name,
                                                    sweep_variables=sweep_vars
                                                ),
                                                EliminationDayInsetAnalyzer(exp_name=expt_name,
                                                                            sweep_variables=sweep_vars
                                                                            ),
                                                EliminationSetDayInsetAnalyzer(exp_name=expt_name,
                                                                               sweep_variables=sweep_vars
                                                                               ),
                                                InsetAnalyzer(exp_name=expt_name,
                                                              sweep_variables=sweep_vars,
                                                              channels=[
                                                                  'Adult Vectors', 'Infectious Vectors',
                                                                  'True Prevalence', 'PfHRP2 Prevalence'
                                                              ]
                                                              )
                                            ]
                                            )
        print(am.experiments)
        am.analyze()
