import numpy as np
import os
import pandas as pd

##
# -------- Set experiments/work items to load

# ---- spatial, classic, GM only, EIR = 10 --> DONE + REDONE
# NOTE THAT THE 6 SWEEP VAR WORK ITEM DOESN'T HAVE ALLELE FREQS (EXP DOES)
# wi_name = 'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne'
# wi_name_sh = 'spatial, classic drive, GM only, EIR = 10'
# wi_names = [
#     # 'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_release_day_release_node_num',
#     'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne',
#     'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newrr0',
#     'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newsne'
# ]
# # AFTER ADDING IN MORE INSET CHART VARS AND RE-RUNNING THE 1ST EXP W/ 4 VARS
# data_dirs = ['Y:\\home\\sleung\\workitems\\dc0\\9b9\\f06\\dc09b9f0-6c25-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\a3e\\2a3\\080\\a3e2a308-031b-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\a93\\a13\\440\\a93a1344-031b-ec11-a9ed-b88303911bc1']
# num_sweep_vars = 4  # 6
# drive_type = 'classic'

# ---- spatial, integral, GM only, EIR = 10 --> DONE + REDONE
# wi_name = 'spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2'
# wi_name_sh = 'spatial, integral drive, GM only, EIR = 10'
# wi_names = ['spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2',
#             'spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newrr20',
#             'spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newse2']
# # AFTER ADDING IN MORE INSET CHART VARS AND RE-RUNNING THE 1ST EXP W/ 4 VARS
# data_dirs = ['Y:\\home\\sleung\\workitems\\fa0\\898\\895\\fa089889-5c25-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\d65\\c35\\7d0\\d65c357d-041b-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\568\\cac\\a30\\568caca3-041b-ec11-a9ed-b88303911bc1']
# num_sweep_vars = 4
# drive_type = 'integral'

# ---- spatial, classic, VC and GM, EIR = 10 --> DONE DONE
# wi_name = 'spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne'
# wi_name_sh = 'spatial, classic drive, VC and GM, EIR = 10'
# wi_names = ['spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne',
#             'spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newrr0',
#             'spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newsne']
# # AFTER ADDING IN MORE INSET CHART VARS AND RE-RUNNING THE 1ST EXP W/ 4 VARS
# data_dirs = ['Y:\\home\\sleung\\workitems\\17c\\74c\\990\\17c74c99-031b-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\1c1\\3aa\\f20\\1c13aaf2-031b-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\33b\\c08\\c30\\33bc08c3-031b-ec11-a9ed-b88303911bc1']
# num_sweep_vars = 4
# drive_type = 'classic'

# ---- spatial, integral, VC and GM, EIR = 10 --> DONE DONE
# wi_name = 'spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2'
# wi_name_sh = 'spatial, integral drive, VC and GM, EIR = 10'
# wi_names = ['spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2',
#             'spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newrr20',
#             'spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newse2']
# # AFTER ADDING IN MORE INSET CHART VARS
# data_dirs = ['Y:\\home\\sleung\\workitems\\24f\\a8c\\e80\\24fa8ce8-041b-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\b6a\\0f7\\0b0\\b6a0f70b-051b-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\e0b\\c2e\\360\\e0bc2e36-051b-ec11-a9ed-b88303911bc1']
# num_sweep_vars = 4
# drive_type = 'integral'

# ---- spatial, classic, GM only, EIR = 30 --> DONE + REDONE
# NOTE THAT THE 6 SWEEP VAR WORK ITEM DOESN'T HAVE ALLELE FREQS (EXP DOES)
# wi_name = 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne'
# wi_name_sh = 'spatial, classic drive, GM only, EIR = 30'
# wi_names = [
#     # 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_release_day_release_node_num',
#     'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne',
#     'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newrr0',
#     'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newsne'
# ]
# # AFTER ADDING IN MORE INSET CHART VARS AND RE-RUNNING THE 1ST EXP W/ 4 VARS
# data_dirs = ['Y:\\home\\sleung\\workitems\\d40\\1fd\\e85\\d401fde8-5a25-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\244\\8a9\\816\\2448a981-6e1a-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\7c4\\563\\fc6\\7c4563fc-6e1a-ec11-a9ed-b88303911bc1']
# num_sweep_vars = 4  # 6
# drive_type = 'classic'

# ---- spatial, integral, GM only, EIR = 30 --> DONE + REDONE
# wi_name = 'spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2'
# wi_name_sh = 'spatial, integral drive, GM only, EIR = 30'
# wi_names = ['spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2',
#             'spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newrr20',
#             'spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newse2']
# # AFTER ADDING IN MORE INSET CHART VARS AND RE-RUNNING THE 1ST EXP W/ 4 VARS
# data_dirs = ['Z:\\home\\sleung\\workitems\\7a1\\d07\\d90\\7a1d07d9-0a1d-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\6c9\\89e\\ae7\\6c989eae-701a-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\380\\e9d\\ac5\\380e9dac-5b25-ec11-9ecd-9440c9bee941']
# num_sweep_vars = 4
# drive_type = 'integral'

# ---- spatial, classic, VC and GM, EIR = 30 --> DONE DONE
# wi_name = 'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne'
# wi_name_sh = 'spatial, classic drive, VC and GM, EIR = 30'
# wi_names = ['spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne',
#             'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newrr0',
#             'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newsne']
# # AFTER ADDING IN MORE INSET CHART VARS
# data_dirs = ['Y:\\home\\sleung\\workitems\\a46\\048\\996\\a4604899-6f1a-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\24c\\b36\\fc6\\24cb36fc-6f1a-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\f67\\e3c\\d56\\f67e3cd5-6f1a-ec11-a9ed-b88303911bc1']
# num_sweep_vars = 4
# drive_type = 'classic'

# ---- spatial, integral, VC and GM, EIR = 30 --> DONE + REDONE
# wi_name = 'spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2'
# wi_name_sh = 'spatial, integral drive, VC and GM, EIR = 30'
# wi_names = ['spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2',
#             'spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newrr20',
#             'spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newse2']
# # AFTER ADDING IN MORE INSET CHART VARS
# data_dirs = ['Y:\\home\\sleung\\workitems\\bc7\\6fd\\0f7\\bc76fd0f-731a-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\94c\\493\\6c7\\94c4936c-731a-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\b81\\67c\\eb5\\b8167ceb-5b25-ec11-9ecd-9440c9bee941']
# num_sweep_vars = 4
# drive_type = 'integral'

# ---- spatial, classic, VC and GM, EIR = 80 --> DONE
# NOTE THAT THE 6 SWEEP VAR WORK ITEM DOESN'T HAVE ALLELE FREQS (EXP DOES)
wi_name = 'spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne'
wi_name_sh = 'spatial, classic drive, VC and GM, EIR = 80'
wi_names = [
    # 'spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_release_day_release_node_num',
    'spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne',
    'spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newrr0',
    'spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newsne'
]
# AFTER ADDING IN MORE INSET CHART VARS AND RE-RUNNING THE 1ST EXP W/ 4 VARS
data_dirs = ['Y:\\home\\sleung\\workitems\\c36\\e72\\a70\\c36e72a7-051b-ec11-a9ed-b88303911bc1',
             'Y:\\home\\sleung\\workitems\\39d\\815\\c40\\39d815c4-061b-ec11-a9ed-b88303911bc1',
             'Y:\\home\\sleung\\workitems\\d0e\\42f\\140\\d0e42f14-061b-ec11-a9ed-b88303911bc1']
num_sweep_vars = 4  # 6
drive_type = 'classic'

# ---- spatial, integral, VC and GM, EIR = 80 --> DONE DONE
# wi_name = 'spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2'
# wi_name_sh = 'spatial, integral drive, VC and GM, EIR = 80'
# wi_names = ['spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2',
#             'spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newrr20',
#             'spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newse2']
# # AFTER ADDING IN MORE INSET CHART VARS
# data_dirs = ['Y:\\home\\sleung\\workitems\\e49\\e74\\0e0\\e49e740e-071b-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\a9e\\e4c\\440\\a9ee4c44-071b-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\2a1\\4de\\7c0\\2a14de7c-071b-ec11-a9ed-b88303911bc1']
# num_sweep_vars = 4
# drive_type = 'integral'

elim_day = 2555  # day on which elim fraction is calculated

file_suffixes = []
for i in range(3):
    file_suffixes.append([])
if (wi_name == 'spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2') \
        or (wi_name == 'spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2') \
        or (wi_name == 'spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2') \
        or (wi_name == 'spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2'):
    # - 1st and 3rd work items
    partition_vars = ['d1']
    partition_vars_vals = [[1, 0.95, 0.9]]
    for partition_vars_val0 in partition_vars_vals[0]:
        fsbegtemp = partition_vars[0] + str(partition_vars_val0)
        file_suffixes[0].append(fsbegtemp)
        file_suffixes[2].append(fsbegtemp)
    # - 2nd work item has no partition vars
elif (wi_name == 'spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2'):
    # - 3rd work item
    partition_vars = ['d1']
    partition_vars_vals = [[1, 0.95, 0.9]]
    for partition_vars_val0 in partition_vars_vals[0]:
        fsbegtemp = partition_vars[0] + str(partition_vars_val0)
        file_suffixes[2].append(fsbegtemp)
    # - 1st and 2nd work items have no partition vars
elif (wi_name == 'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne') \
        or (wi_name == 'spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne') \
        or (wi_name == 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne') \
        or (wi_name == 'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne') \
        or (wi_name == 'spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne'):
    # - 1st, 2nd, and 3rd work items have no partition vars
    pass

# -------- Load data
dfi = pd.DataFrame()
dfa = pd.DataFrame()
dfe = pd.DataFrame()
dfed = pd.DataFrame()
for iwn in range(len(file_suffixes)):
    file_suffixesnow = file_suffixes[iwn]
    data_dirnow = data_dirs[iwn]
    wi_namenow = wi_names[iwn]
    if len(file_suffixesnow) > 0:
        for file_suffix in file_suffixesnow:
            filei = os.path.join(data_dirnow, wi_namenow + '_inset_data_' + file_suffix + '.csv')
            filea = os.path.join(data_dirnow, wi_namenow + '_spatial_avg_allele_freqs_' + file_suffix + '.csv')
            filee = os.path.join(data_dirnow, wi_namenow + '_inset_data_elim_day_'
                                 + str(elim_day) + '_indiv_sims_' + file_suffix + '.csv')
            fileed = os.path.join(data_dirnow,
                                  wi_namenow + '_inset_data_elim_day_number_indiv_sims_' + file_suffix + '.csv')
            dfi = dfi.append(pd.read_csv(filei))
            dfa = dfa.append(pd.read_csv(filea))
            dfe = dfe.append(pd.read_csv(filee))
            dfed = dfed.append(pd.read_csv(fileed))
    else:
        dfi = dfi.append(pd.read_csv(os.path.join(data_dirnow, wi_namenow + '_inset_data.csv')))
        dfa = dfa.append(pd.read_csv(os.path.join(data_dirnow, wi_namenow + '_spatial_avg_allele_freqs.csv')))
        dfe = dfe.append(pd.read_csv(
            os.path.join(data_dirnow, wi_namenow + '_inset_data_elim_day_' + str(elim_day) + '_indiv_sims.csv')))
        dfed = dfed.append(
            pd.read_csv(os.path.join(data_dirnow, wi_namenow + '_inset_data_elim_day_number_indiv_sims.csv')))

# -------- Clean up data
dfe = dfe.drop(columns=['Daily_EIR_elim', 'New_Clinical_Cases_elim', 'Run_Number'])
dfed = dfed.drop(columns=['Run_Number'])
if drive_type == 'integral':
    dfi = dfi[dfi['rr20'] != 0.2]
    dfa = dfa[dfa['rr20'] != 0.2]
    dfe = dfe[dfe['rr20'] != 0.2]
    dfed = dfed[dfed['rr20'] != 0.2]
    if num_sweep_vars == 4:
        dfi = dfi[['Time', 'rc', 'd1', 'rr20', 'se2',
                   'PfHRP2 Prevalence', 'PfHRP2 Prevalence_std',
                   'True Prevalence', 'True Prevalence_std',
                   'Adult Vectors', 'Adult Vectors_std',
                   'Infectious Vectors', 'Infectious Vectors_std',
                   # 'Daily EIR', 'Daily EIR_std'
                   ]]
elif drive_type == 'classic':
    dfi = dfi[dfi['rr0'] != 0.2]
    dfa = dfa[dfa['rr0'] != 0.2]
    dfe = dfe[dfe['rr0'] != 0.2]
    dfed = dfed[dfed['rr0'] != 0.2]
    if num_sweep_vars == 6:
        dfi.rename(columns={'release_day': 'rd', 'num_nodes': 'nn'}, inplace=True)
        dfa.rename(columns={'release_day': 'rd', 'num_nodes': 'nn'}, inplace=True)
        dfe.rename(columns={'release_day': 'rd', 'num_nodes': 'nn'}, inplace=True)
        dfed.rename(columns={'release_day': 'rd', 'num_nodes': 'nn'}, inplace=True)
        dfi = dfi[['Time', 'rc', 'd', 'rr0', 'sne', 'rd', 'nn',
                   'PfHRP2 Prevalence', 'PfHRP2 Prevalence_std',
                   'True Prevalence', 'True Prevalence_std',
                   'Adult Vectors', 'Adult Vectors_std',
                   'Infectious Vectors', 'Infectious Vectors_std',
                   # 'Daily EIR', 'Daily EIR_std'
                   ]]
    elif num_sweep_vars == 4:
        dfi = dfi[['Time', 'rc', 'd', 'rr0', 'sne',
                   'PfHRP2 Prevalence', 'PfHRP2 Prevalence_std',
                   'True Prevalence', 'True Prevalence_std',
                   'Adult Vectors', 'Adult Vectors_std',
                   'Infectious Vectors', 'Infectious Vectors_std',
                   # 'Daily EIR', 'Daily EIR_std'
                   ]]

# -------- Save out data
dfi.to_csv('csvs/dfi_' + wi_name + '.csv', index=False)
dfa.to_csv('csvs/dfa_' + wi_name + '.csv', index=False)
dfe.to_csv('csvs/dfe_' + wi_name + '.csv', index=False)
dfed.to_csv('csvs/dfed_' + wi_name + '.csv', index=False)
