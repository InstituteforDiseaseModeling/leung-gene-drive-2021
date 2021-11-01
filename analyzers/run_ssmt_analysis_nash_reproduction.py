from simtools.Managers.WorkItemManager import WorkItemManager
from simtools.SetupParser import SetupParser
from simtools.AssetManager.FileList import FileList

# wi_name = "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne"
# wi_name = "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne"
# wi_name = "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne"
# wi_name = "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newrr0"
# wi_name = "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newrr20"
# wi_name = "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newrr20"
# wi_name = "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newrr20"
# wi_name = "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newrr0"
# wi_name = "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newrr0"
# wi_name = "spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newsne"
# wi_name = "spatialinside_classic3allele_GM_only_aEIR80_sweep_rc_d_rr0_sne_newsne"
# wi_name = "spatialinside_integral2l4a_GM_only_aEIR80_sweep_rc_d1_rr20_se2_newse2"
# wi_name = "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newse2"
# wi_name = "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newse2"
# wi_name = "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newse2"
# wi_name = "spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newsne"
# wi_name = "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newsne"
# wi_name = "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newsne"
# wi_name = "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newsne"
# wi_name = "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newse2"
# wi_name = "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newse2"
# wi_name = "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newrr20"
# wi_name = "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newrr20"
# wi_name = "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newrr0"
# wi_name = "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newrr0"
# wi_name = "spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_sne_newsne"
# wi_name = "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_sne_newsne"
# wi_name = "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_se2_newse2"
# wi_name = "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_se2_newse2"
# wi_name = "spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne"
# wi_name = "spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne"
wi_name = "spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2"
# wi_name = "spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2"
# wi_name = "spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2"
# wi_name = "spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2"
# wi_name = "spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2"
command = "python3 run_analysis_nash_reproduction_spatial.py"

# ------------------------------------
user_files = FileList(root='analyzers_large_data')

if __name__ == "__main__":
    SetupParser.default_block = "HPC"
    SetupParser.init()

    wim = WorkItemManager(item_name=wi_name, command=command, user_files=user_files,
                          docker_image="docker-production.packages.idmod.org/dtk-tools/dtk-tools-feather-tqdm-hdf5:1.0.1",
                          related_experiments=["2810b14b-b125-ec11-9ecd-9440c9bee941"])
    wim.comps_env = 'Calculon'
    wim.execute(True)
