
# COVERAGE ----------------------------------------------------------------


for(nobs in c(10,100,1000)){
  umm = readRDS(paste0("./analysis/supp_material/simulation_output/coverage_PI_61_",nobs,".rds"))
  print(apply(umm, 2, mean)) 
}

for(nobs in c(10,100,1000)){
  umm = readRDS(paste0("./analysis/supp_material/simulation_output/coverage_EI_61_",nobs,".rds"))
  print(apply(umm, 2, mean)) 
}

for(nobs in c(10,100,1000)){
  umm = readRDS(paste0("./analysis/supp_material/simulation_output/coverage_UCB_61_",nobs,".rds"))
  print(apply(umm, 2, mean)) 
}



# OMITTED -----------------------------------------------------------------

for(nobs in c(10,100,1000)){
  umm = readRDS(paste0("./analysis/supp_material/simulation_output/PI_EI_omitted_1_",nobs,".rds"))
  print(apply(umm, 2, mean, na.rm = T)) 
}

