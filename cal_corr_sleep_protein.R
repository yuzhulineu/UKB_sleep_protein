### cal the association between sleep traits and proteins

cl<-makeCluster(15)
registerDoParallel(cl)
result_sleep_protein<-list()

for (kk in 1:8) {
  result_value_protein<-data.frame(matrix(NA,ncol=4,nrow=2923))
  sleep_traits1<-sleep_final[,kk]
  result_protein<-foreach(
    ii=1:2922,.combine = rbind) %dopar% calcorr_sleep_protein(ii,sleep_traits1,result_value_protein)
  result_sleep_protein[[kk]]<-result_protein
  print(kk)
}


calcorr_sleep_protein<-function(kk,sleep_traits1,result_value_protein) {
  protein_final<-sleep_protein_cov_batch[,kk+1]

  fit1<-lm(sleep_traits1~protein_final+Age+Sex+Townsand+BMI+Qualification+Smoking+Alcohal
           +Batch+gap_day,data=sleep_protein_cov_batch)
  result_summary<-summary(fit1)
  result_value_protein[kk,]<-result_summary$coefficients[2,]
}
