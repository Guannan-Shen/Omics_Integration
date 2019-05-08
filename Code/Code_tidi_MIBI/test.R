nb_fam_diagnosis_age_sex <- all_fam %>%      ## tidy_MIBI set
  NB_mods(rank = "Family", ## Rank of taxa we want to model
          diagnosis, Age, Sex       ## The covariates in our model (Group + EWG)
  )


# vars = c("CD4_cnt")
# 
# nb_fam <- all_fam %>%      ## tidy_MIBI set
#   NB_mods(rank = "Family", ## Rank of taxa we want to model
#           eval(parse(text = "CD4_Cnt")), Age, Sex       ## The covariates in our model (Group + EWG)
#   )
# 
# eval(as.name("data1"))
# 
# eval()
# 
# Age
getwd()
save(nb_fam_diagnosis_age_sex, 
     file = "C:/Users/hithr/Documents/Stats/gitlab/tidi_MIBI/DataProcessed/nb_model/nb_fam_diagnosis_age_sex.RData" )