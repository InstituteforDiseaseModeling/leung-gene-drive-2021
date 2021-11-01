import os
import pandas as pd

##
# -----------------------------------------
# Set experiments/work items to load
# -----------------------------------------

# ---- spatial, classic, GM only, EIR = 10
# wi_name = 'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne'
# wi_name_sh = 'spatial, classic drive, GM only, EIR = 10'
# wi_names = [
#     'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne',
#     'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newrr0',
#     'spatialinside_classic3allele_GM_only_aEIR10_sweep_rc_d_rr0_sne_newsne'
# ]
# data_dirs = ['Y:\\home\\sleung\\workitems\\dc0\\9b9\\f06\\dc09b9f0-6c25-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\091\\bd2\\ac1\\091bd2ac-1826-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\1c6\\834\\260\\1c683426-0e26-ec11-9ecd-9440c9bee941']
# num_sweep_vars = 4
# drive_type = 'classic'

# ---- spatial, integral, GM only, EIR = 10
# wi_name = 'spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2'
# wi_name_sh = 'spatial, integral drive, GM only, EIR = 10'
# wi_names = ['spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2',
#             'spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newrr20',
#             'spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2_newse2']
# data_dirs = ['Y:\\home\\sleung\\workitems\\fa0\\898\\895\\fa089889-5c25-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\d65\\c35\\7d0\\d65c357d-041b-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\568\\cac\\a30\\568caca3-041b-ec11-a9ed-b88303911bc1']
# num_sweep_vars = 4
# drive_type = 'integral'

# ---- spatial, classic, VC and GM, EIR = 10
# wi_name = 'spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne'
# wi_name_sh = 'spatial, classic drive, VC and GM, EIR = 10'
# wi_names = ['spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne',
#             'spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newrr0',
#             'spatialinside_classic3allele_VC_and_GM_aEIR10_sweep_rc_d_rr0_sne_newsne']
# data_dirs = ['Y:\\home\\sleung\\workitems\\044\\a98\\9c0\\044a989c-0a27-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\d3a\\6ca\\5b3\\d3a6ca5b-3526-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\801\\54e\\0d3\\80154e0d-3526-ec11-9ecd-9440c9bee941']
# num_sweep_vars = 4
# drive_type = 'classic'

# ---- spatial, integral, VC and GM, EIR = 10
# wi_name = 'spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2'
# wi_name_sh = 'spatial, integral drive, VC and GM, EIR = 10'
# wi_names = ['spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2',
#             'spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newrr20',
#             'spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2_newse2']
# data_dirs = ['Y:\\home\\sleung\\workitems\\793\\1c6\\f37\\7931c6f3-7426-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\b6a\\0f7\\0b0\\b6a0f70b-051b-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\e0b\\c2e\\360\\e0bc2e36-051b-ec11-a9ed-b88303911bc1']
# num_sweep_vars = 4
# drive_type = 'integral'

# ---- spatial, classic, GM only, EIR = 30
# wi_name = 'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne'
# wi_name_sh = 'spatial, classic drive, GM only, EIR = 30'
# wi_names = [
#     'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne',
#     'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newrr0',
#     'spatialinside_classic3allele_GM_only_aEIR30_sweep_rc_d_rr0_sne_newsne'
# ]
# data_dirs = ['Y:\\home\\sleung\\workitems\\d40\\1fd\\e85\\d401fde8-5a25-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\155\\f84\\457\\155f8445-7126-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\9b7\\3b3\\9c0\\9b73b39c-0c26-ec11-9ecd-9440c9bee941']
# num_sweep_vars = 4
# drive_type = 'classic'

# ---- spatial, integral, GM only, EIR = 30
# wi_name = 'spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2'
# wi_name_sh = 'spatial, integral drive, GM only, EIR = 30'
# wi_names = ['spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2',
#             'spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newrr20',
#             'spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2_newse2']
# data_dirs = ['Y:\\home\\sleung\\workitems\\c5c\\376\\3d0\\c5c3763d-0d26-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\6c9\\89e\\ae7\\6c989eae-701a-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\380\\e9d\\ac5\\380e9dac-5b25-ec11-9ecd-9440c9bee941']
# num_sweep_vars = 4
# drive_type = 'integral'

# ---- spatial, classic, VC and GM, EIR = 30
# wi_name = 'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne'
# wi_name_sh = 'spatial, classic drive, VC and GM, EIR = 30'
# wi_names = ['spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne',
#             'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newrr0',
#             'spatialinside_classic3allele_VC_and_GM_aEIR30_sweep_rc_d_rr0_sne_newsne']
# data_dirs = ['Y:\\home\\sleung\\workitems\\efc\\615\\fd7\\efc615fd-7326-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\50d\\1d7\\8c7\\50d1d78c-7126-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\89c\\c04\\a20\\89cc04a2-0d26-ec11-9ecd-9440c9bee941']
# num_sweep_vars = 4
# drive_type = 'classic'

# ---- spatial, integral, VC and GM, EIR = 30
# wi_name = 'spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2'
# wi_name_sh = 'spatial, integral drive, VC and GM, EIR = 30'
# wi_names = ['spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2',
#             'spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newrr20',
#             'spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2_newse2']
# data_dirs = ['Y:\\home\\sleung\\workitems\\3f0\\f8e\\247\\3f0f8e24-7226-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\94c\\493\\6c7\\94c4936c-731a-ec11-a9ed-b88303911bc1',
#              'Y:\\home\\sleung\\workitems\\b81\\67c\\eb5\\b8167ceb-5b25-ec11-9ecd-9440c9bee941']
# num_sweep_vars = 4
# drive_type = 'integral'

# ---- spatial, classic, VC and GM, EIR = 80
# wi_name = 'spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne'
# wi_name_sh = 'spatial, classic drive, VC and GM, EIR = 80'
# wi_names = [
#     'spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne',
#     'spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newrr0',
#     'spatialinside_classic3allele_VC_and_GM_aEIR80_sweep_rc_d_rr0_sne_newsne'
# ]
# data_dirs = ['Y:\\home\\sleung\\workitems\\8b6\\f4d\\360\\8b6f4d36-0c26-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\18d\\979\\6e0\\18d9796e-0b27-ec11-9ecd-9440c9bee941',
#              'Y:\\home\\sleung\\workitems\\9bd\\6c2\\2d0\\9bd6c22d-0b27-ec11-9ecd-9440c9bee941']
# num_sweep_vars = 4
# drive_type = 'classic'

# ---- spatial, integral, VC and GM, EIR = 80
wi_name = 'spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2'
wi_name_sh = 'spatial, integral drive, VC and GM, EIR = 80'
wi_names = ['spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2',
            'spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newrr20',
            'spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2_newse2']
data_dirs = ['Y:\\home\\sleung\\workitems\\7f8\\630\\fe0\\7f8630fe-0b27-ec11-9ecd-9440c9bee941',
             'Y:\\home\\sleung\\workitems\\a9e\\e4c\\440\\a9ee4c44-071b-ec11-a9ed-b88303911bc1',
             'Y:\\home\\sleung\\workitems\\2a1\\4de\\7c0\\2a14de7c-071b-ec11-a9ed-b88303911bc1']
num_sweep_vars = 4
drive_type = 'integral'

# -----------------------------------------
# Gather all csvs and save out one
# -----------------------------------------
elim_day = 2555  # day on which elim fraction is calculated

file_suffixes = []
for i in range(3):
    file_suffixes.append([])
if (wi_name == 'spatialinside_integral2l4a_VC_and_GM_aEIR80_sweep_rc_d1_rr20_se2') \
        or (wi_name == 'spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2') \
        or (wi_name == 'spatialinside_integral2l4a_GM_only_aEIR30_sweep_rc_d1_rr20_se2') \
        or (wi_name == 'spatialinside_integral2l4a_VC_and_GM_aEIR10_sweep_rc_d1_rr20_se2') \
        or (wi_name == 'spatialinside_integral2l4a_VC_and_GM_aEIR30_sweep_rc_d1_rr20_se2'):
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
dfi.to_csv('../csvs/dfi_' + wi_name + '.csv', index=False)
dfa.to_csv('../csvs/dfa_' + wi_name + '.csv', index=False)
dfe.to_csv('../csvs/dfe_' + wi_name + '.csv', index=False)
dfed.to_csv('../csvs/dfed_' + wi_name + '.csv', index=False)
