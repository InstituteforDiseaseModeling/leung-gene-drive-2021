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
    #     "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_release_node_num":
    #         "0d04ced1-78f1-eb11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_release_day_release_node_num":
    #         "60cdc1b7-79f1-eb11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_release_day_release_node_num":
    #         "ee6880ac-4dee-eb11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2":
    #         "7fdeaca5-eefe-eb11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2":
    #         "187be08e-f501-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2":
    #         "655f1942-5004-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2":
    #         "c5865cba-5004-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2":
    #         "459ec043-5104-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne":
    #         "f7ad7e26-a404-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne":
    #         "8446bdd2-a404-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_release_number":
    #         "686080e8-8b0a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_release_number":
    #         "efd2be68-8b0a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_release_number":
    #         "b8b506c3-8d0a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR30_sweep_release_number":
    #         "30f67e66-8d0a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_se2_newse2":
    #         "0194b407-e60a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_se2_newse2":
    #         "f90e04ab-e60a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_sne_newsne":
    #         "81a0d1c8-c30a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_sne_newsne":
    #         "c7809967-c40a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newrr0":
    #         "c9d8e922-860a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newrr0":
    #         "5e309994-860a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newrr20":
    #         "32d788f3-b70a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newrr20":
    #         "b2bdccfc-b80a-ec11-a9ed-b88303911bc1"
    # }

    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newse2":
    #         "7042c63e-0f12-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newse2":
    #         "82b2b8f0-1012-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newsne":
    #         "42762228-1212-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newsne":
    #         "04d92676-1212-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newsne":
    #         "89885530-6812-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newsne":
    #         "20c068d2-6812-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newse2":
    #         "05e863ae-7f12-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newse2":
    #         "f5ab2c5b-8012-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newse2":
    #         "b4b8bbb4-8012-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR80_sweep_rc_d1_rr20_se2_newse2":
    #         "79961938-8112-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR80_sweep_rc_d_rr0_sne_newsne":
    #         "171581cb-8112-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newsne":
    #         "d369c907-8212-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newrr0":
    #         "65c56359-8c12-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newrr0":
    #         "cbd3fab9-8c12-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newrr20":
    #         "81a418c1-a212-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newrr20":
    #         "4578e934-a312-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newrr20":
    #         "e91d8c88-a312-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newrr0":
    #         "92711bee-a312-ec11-a9ed-b88303911bc1"
    # }

    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newrr0":
    #         "5e309994-860a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newsne":
    #         "04d92676-1212-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne":
    #         "f7ad7e26-a404-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newsne":
    #         "42762228-1212-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newrr0":
    #         "c9d8e922-860a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2":
    #         "7fdeaca5-eefe-eb11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newrr20":
    #         "b2bdccfc-b80a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newse2":
    #         "82b2b8f0-1012-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2":
    #         "655f1942-5004-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newrr20":
    #         "32d788f3-b70a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newse2":
    #         "7042c63e-0f12-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne":
    #         "1a250220-a71a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne":
    #         "292c08d2-a61a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newrr0":
    #         "cbd3fab9-8c12-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newsne":
    #         "20c068d2-6812-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne":
    #         "8446bdd2-a404-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newsne":
    #         "89885530-6812-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newrr0":
    #         "65c56359-8c12-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2":
    #         "187be08e-f501-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newrr20":
    #         "4578e934-a312-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newse2":
    #         "05e863ae-7f12-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2":
    #         "c5865cba-5004-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newrr20":
    #         "81a418c1-a212-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newse2":
    #         "f5ab2c5b-8012-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne":
    #         "b17c24d9-a71a-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newsne":
    #         "d369c907-8212-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newrr0":
    #         "92711bee-a312-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2":
    #         "459ec043-5104-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newrr20":
    #         "e91d8c88-a312-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newse2":
    #         "b4b8bbb4-8012-ec11-a9ed-b88303911bc1"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rr0_sne_rc0.9_d0.95_50seeds":
    #         "47bea906-ec1f-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rr0_sne_rc0.8_d1_50seeds":
    #         "7040e599-eb1f-ec11-9ecd-9440c9bee941"
    # }

    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne":
    #         "f5a059df-ae21-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newse2":
    #         "9250cec6-af21-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newse2":
    #         "79f3816b-b021-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2":
    #         "87952b20-b221-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne":
    #         "a661f6b1-5525-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne":
    #         "0ee573f0-5625-ec11-9ecd-9440c9bee941"
    # }

    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newsne":
    #         "e265629a-a725-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2":
    #         "51b0ed6e-a825-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newsne":
    #         "e6a98024-a925-ec11-9ecd-9440c9bee941"
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
    #     "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newsne":
    #         "1b495c58-ac25-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newrr0":
    #         "c004e0ed-ac25-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newrr0":
    #         "253b91f5-2826-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newrr0":
    #         "755c1724-3626-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2":
    #         "1e2dd6d7-2926-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne":
    #         "4abc9dfa-2a26-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2":
    #         "1a3f69cf-b025-ec11-9ecd-9440c9bee941"
    # }
    # experiments = {
    #     "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne":
    #         "940c6484-3626-ec11-9ecd-9440c9bee941"
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
