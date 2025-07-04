setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/Datasets/StatsSampleNames/')

prefix_namerica = 'minCount2.namerica_ref_'
prefix_europe = 'minCount2.europe_ref_'

array_chroms = c('LG14_SUPER_15',
                 'LG9_SUPER_8',
                 'LG13_SUPER_17',
                 'LG5_SUPER_9',
                 'LG6_SUPER_11',
                 'LG11_SUPER_16',
                 'LG16_SUPER_12',
                 'LG15_SUPER_13',
                 'LG4_SUPER_7',
                 'LG10_SUPER_14',
                 'LG7_SUPER_10',
                 'LG17_SUPER_6',
                 'LG8_SUPER_5',
                 'LG12_SUPER_3',
                 'LG3_SUPER_4',
                 'LG1_SUPER_1',
                 'LG2_SUPER_2')

array_namerica_samples = c('9NSC',
                           '8NJY',
                           '7NFD',
                           '14SJS',
                           '13SHO',
                           '10NUK',
                           'PTY_BLMP2',
                           'PLIT_BLMP3',
                           'PLIB_BLMP4',
                           'PHOLY_BLMP5',
                           'PFOX_BLMP6',
                           'PBY_BLMP1',
                           'PLY-high',
                           'RED',
                           'OAK',
                           'SM-BOO-high-merged',
                           'SM-BOO-low-merged')

array_europe_samples = c('3GAL',
                         '2FAR',
                         '1BER',
                         '18VAD',
                         '17TFA',
                         '16STK',
                         '15SKA',
                         '12POR',
                         'VA-3191-TSWANG-low_S24_L003',
                         'VA-3191-TSWANG-high_S25_L003',
                         'VA-3191-SPCAN-low_S22_L003',
                         'VA-3191-SPCAN-high_S23_L003',
                         'VA-3191-Bodo-low_S26_L003',
                         'VA-3191-Bodo-high_S27_L003',
                         'EAS-high',
                         'EAS-low',
                         'T_NC_NO-high',
                         'VEN',
                         'SPn_Crab_High',
                         'SPn_Wave_Low',
                         'Fr_Crab_Low',
                         'Fr_Wave_High',
                         'SPs_Crab_High',
                         'SPs_Wave_Low',
                         'Uke_Wave_High',
                         'Uke_Crab_Low',
                         'SWs_Wave',
                         'SWs_Crab',
                         'Ukw_Wave_High',
                         'Ukw_Crab_Low',
                         'SM-SWn3_Wave-merged',
                         'SM-SWn3_Crab-merged',
                         'SM-SWn5_Crab-merged',
                         'SM-SWn5_Wave-merged')



for ( chrom in array_chroms) {
  tmp_array_chrom_names_namerica = paste0(prefix_namerica, chrom, '_java.', seq(1:length(array_namerica_samples)))
  tmp_array_chrom_names_europe = paste0(prefix_europe, chrom, '_java.', seq(1:length(array_europe_samples)))
  
  df_names_namerica = data.frame(tmp_array_chrom_names_namerica, array_namerica_samples)
  df_names_europe = data.frame(tmp_array_chrom_names_europe, array_europe_samples)
  
  file_name_namerica = paste('Namerica', chrom,'rename.tsv', sep = '.')
  file_name_europe = paste('Europe', chrom,'rename.tsv', sep = '.')
  
  print(file_name_namerica)
  
  write.table(df_names_namerica, file = file_name_namerica, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  write.table(df_names_europe, file = file_name_europe, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}
