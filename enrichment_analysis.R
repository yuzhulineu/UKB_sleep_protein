# perform the enrichment analysis
library(gprofiler2)

####################short
p_value_shortsleep<-p_value_proteinfinal1$short.sleep
#index_shortsleep<-p_value_shortsleep<0.01/1463
index_shortsleep<-which(p_value_shortsleep<0.05/(2922*8))
index_negative<-which(beta_value_proteinfinal1$short.sleep<0)
index_negativeprotein<-intersect(index_shortsleep,index_negative)
protein_short_n<-beta_value_proteinfinal1$Assay[index_negativeprotein]


enrich_result_short1 <- gost(query = protein_short_n,
                             organism = 'hsapiens', ordered_query = FALSE,
                             domain_scope = "custom",custom_bg = beta_value_proteinfinal1$Assay,
                             multi_query = FALSE, exclude_iea = FALSE,
                             measure_underrepresentation = FALSE, evcodes = FALSE,
                             sources = c('GO','KEGG','REAC','WP'),
                             significant = TRUE, user_threshold = 0.01,correction_method = "fdr")#correction_method = "g_SCS"

#gostplot(enrich_result_short1, interactive = TRUE)
enrich_result_short_negative<-enrich_result_short1$result

######################long sleep
p_value_longsleep<-p_value_proteinfinal1$long.sleep
#index_shortsleep<-p_value_shortsleep<0.01/1463
index_longsleep<-which(p_value_longsleep<0.05/(2922*8))
index_positive<-which(beta_value_proteinfinal1$long.sleep>0)
index_positiveprotein<-intersect(index_longsleep,index_positive)
protein_long_n<-beta_value_proteinfinal1$Assay[index_positiveprotein]


enrich_result_long1 <- gost(query = protein_long_n,
                             organism = 'hsapiens', ordered_query = FALSE,
                             domain_scope = "custom",custom_bg = beta_value_proteinfinal1$Assay,
                             multi_query = FALSE, exclude_iea = FALSE,
                             measure_underrepresentation = FALSE, evcodes = FALSE,
                            sources = c('GO','KEGG','REAC','WP'),
                             significant = TRUE, user_threshold = 0.01,correction_method = "fdr")#correction_method = "g_SCS"

#gostplot(enrich_result_long1, interactive = TRUE)
enrich_result_long_positive<-enrich_result_long1$result

#########long sleep negative
p_value_longsleep<-p_value_proteinfinal1$long.sleep
#index_shortsleep<-p_value_shortsleep<0.01/1463
index_longsleep<-which(p_value_longsleep<0.05/(2922*8))
index_negative<-which(beta_value_proteinfinal1$long.sleep<0)
index_negativeprotein<-intersect(index_longsleep,index_negative)
protein_long_n<-beta_value_proteinfinal1$Assay[index_negativeprotein]


enrich_result_long1 <- gost(query = protein_long_n,
                            organism = 'hsapiens', ordered_query = FALSE,
                            domain_scope = "custom",custom_bg = beta_value_proteinfinal1$Assay,
                            multi_query = FALSE, exclude_iea = FALSE,
                            measure_underrepresentation = FALSE, evcodes = FALSE,
                            sources = c('GO','KEGG','REAC','WP'),
                            significant = TRUE, user_threshold = 0.01,correction_method = "fdr")#correction_method = "g_SCS"

#gostplot(enrich_result_long1, interactive = TRUE)
enrich_result_long_negative<-enrich_result_long1$result

######################nap during day
p_value_nap<-p_value_proteinfinal1$nap_during_day
#index_shortsleep<-p_value_shortsleep<0.01/1463
index_nap<-which(p_value_nap<0.05/(2922*8))
index_positive<-which(beta_value_proteinfinal1$nap_during_day>0)
index_positiveprotein<-intersect(index_nap,index_positive)
protein_nap_n<-beta_value_proteinfinal1$Assay[index_positiveprotein]


enrich_result_nap1 <- gost(query = protein_nap_n,
                           organism = 'hsapiens', ordered_query = FALSE,
                           domain_scope = "custom",custom_bg = beta_value_proteinfinal1$Assay,
                           multi_query = FALSE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = FALSE,
                           sources = c('GO','KEGG','REAC','WP'),
                           significant = TRUE, user_threshold = 0.01,correction_method = "fdr")#correction_method = "g_SCS"

#gostplot(enrich_result_nap1, interactive = TRUE)
enrich_result_nap_positive<-enrich_result_nap1$result

############ nap during day negative
p_value_nap<-p_value_proteinfinal1$nap_during_day
#index_shortsleep<-p_value_shortsleep<0.01/1463
index_nap<-which(p_value_nap<0.05/(2922*8))
index_negative<-which(beta_value_proteinfinal1$nap_during_day<0)
index_negativeprotein<-intersect(index_nap,index_negative)
protein_nap_n<-beta_value_proteinfinal1$Assay[index_negativeprotein]


enrich_result_nap1 <- gost(query = protein_nap_n,
                           organism = 'hsapiens', ordered_query = FALSE,
                           domain_scope = "custom",custom_bg = beta_value_proteinfinal1$Assay,
                           multi_query = FALSE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = FALSE,
                           sources = c('GO','KEGG','REAC','WP'),
                           significant = TRUE, user_threshold = 0.01,correction_method = "fdr")#correction_method = "g_SCS"

#gostplot(enrich_result_nap1, interactive = TRUE)
enrich_result_nap_negative<-enrich_result_nap1$result

###################get up in morning
p_value_getup<-p_value_proteinfinal1$getting_up_in_morning
#index_getup<-p_value_getup<0.01/1463
index_getup<-which(p_value_getup<0.05/(2922*8))
index_negative<-which(beta_value_proteinfinal1$getting_up_in_morning<0)
index_negativeprotein<-intersect(index_getup,index_negative)
protein_getup_n<-beta_value_proteinfinal1$Assay[index_negativeprotein]


enrich_result_getup1 <- gost(query = protein_getup_n,
                             organism = 'hsapiens', ordered_query = FALSE,
                             domain_scope = "custom",custom_bg = beta_value_proteinfinal1$Assay,
                             multi_query = FALSE, exclude_iea = FALSE,
                             measure_underrepresentation = FALSE, evcodes = FALSE,
                             sources = c('GO','KEGG','REAC','WP'),
                             significant = TRUE, user_threshold = 0.01,correction_method = "fdr")#correction_method = "g_SCS"

#gostplot(enrich_result_getup1, interactive = TRUE)
enrich_result_getup_negative<-enrich_result_getup1$result

####################positive
index_positive<-which(beta_value_proteinfinal1$getting_up_in_morning>0)
index_positiveprotein<-intersect(index_getup,index_positive)
protein_getup_n<-beta_value_proteinfinal1$Assay[index_positiveprotein]


enrich_result_getup1 <- gost(query = protein_getup_n,
                             organism = 'hsapiens', ordered_query = FALSE,
                             domain_scope = "custom",custom_bg = beta_value_proteinfinal1$Assay,
                             multi_query = FALSE, exclude_iea = FALSE,
                             measure_underrepresentation = FALSE, evcodes = FALSE,
                             sources = c('GO','KEGG','REAC','WP'),
                             significant = TRUE, user_threshold = 0.01,correction_method = "fdr")#correction_method = "g_SCS"

#gostplot(enrich_result_getup1, interactive = TRUE)
enrich_result_getup_positive<-enrich_result_getup1$result

######################daytime dozing
p_value_daytimedozing<-p_value_proteinfinal1$daytime_dozing
#index_shortsleep<-p_value_shortsleep<0.01/1463
index_daytimedozing<-which(p_value_daytimedozing<0.05/(2922*8))
index_positive<-which(beta_value_proteinfinal1$daytime_dozing>0)
index_positiveprotein<-intersect(index_daytimedozing,index_positive)
protein_daytimedozing_n<-beta_value_proteinfinal1$Assay[index_positiveprotein]


enrich_result_daytimedozing1 <- gost(query = protein_daytimedozing_n,
                                     organism = 'hsapiens', ordered_query = FALSE,
                                     domain_scope = "custom",custom_bg = beta_value_proteinfinal1$Assay,
                                     multi_query = FALSE, exclude_iea = FALSE,
                                     measure_underrepresentation = FALSE, evcodes = FALSE,
                                     sources = c('GO','KEGG','REAC','WP'),
                                     significant = TRUE, user_threshold = 0.01,correction_method = "fdr")#correction_method = "g_SCS"

#gostplot(enrich_result_daytimedozing1, interactive = TRUE)
enrich_result_daytimedozing_positive<-enrich_result_daytimedozing1$result

############################### negative
p_value_daytimedozing<-p_value_proteinfinal1$daytime_dozing
#index_shortsleep<-p_value_shortsleep<0.01/1463
index_daytimedozing<-which(p_value_daytimedozing<0.05/(2922*8))
index_negative<-which(beta_value_proteinfinal1$daytime_dozing<0)
index_negativeprotein<-intersect(index_daytimedozing,index_negative)
protein_daytimedozing_n<-beta_value_proteinfinal1$Assay[index_negativeprotein]


enrich_result_daytimedozing1 <- gost(query = protein_daytimedozing_n,
                                     organism = 'hsapiens', ordered_query = FALSE,
                                     domain_scope = "custom",custom_bg = beta_value_proteinfinal1$Assay,
                                     multi_query = FALSE, exclude_iea = FALSE,
                                     measure_underrepresentation = FALSE, evcodes = FALSE,
                                     sources = c('GO','KEGG','REAC','WP'),
                                     significant = TRUE, user_threshold = 0.01,correction_method = "fdr")#correction_method = "g_SCS"

#gostplot(enrich_result_daytimedozing1, interactive = TRUE)
enrich_result_daytimedozing_negative<-enrich_result_daytimedozing1$result
