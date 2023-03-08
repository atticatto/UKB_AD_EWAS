#Time-to-event analysis
library(data.table)
library(dplyr)

library(survival)
library(ggfortify)
library(survminer)
library(forcats)

#####0.Prepare survival data #####
base_surv  <- read.csv('/AD_phewas/baseline_table.csv')

######0.1 read covar table######
phesant_covar <- fread('/AD_phewas/phesant_covar.csv')
base_surv <- left_join(base_surv,phesant_covar,by='eid')
base_surv$edu_cat <- factor(base_surv$edu_cat)
str(base_surv$edu_cat)

######0.2 sample filter######
analyzetable <- base_surv
attach(analyzetable)
analyzetable$Age_at_FUend <-Age_at_baseline+AD_day/365 
detach(analyzetable)
#筛随访时限
table(analyzetable$Age_at_baseline)
analyzetable<- filter(analyzetable,!is.na(AD_day),AD_year!=2037,(AD_status==1|AD_status==2&Age_at_FUend>50))
median(analyzetable$AD_day)
quantile(analyzetable$AD_day)
table(analyzetable$Age_at_baseline)
hist(analyzetable$AD_day,breaks = 100)
analyzetable <- analyzetable[,-which(colnames(analyzetable)=='Age_at_FUend')]
analyzetable_backup <- analyzetable#backup

######0.3 read variable table ######
temp <- read.csv('/AD phewas/baseline_table_colnames.csv')

#============================================================================#
#============================================================================#

#####COX analysis#####
#NOTE# Populatoin：Full sample，Follow-up>=3 years
#NOTE# Model：FML <- as.formula(paste0('Surv(AD_day,AD_status)~',paste(i,'Age_at_baseline','sex’，‘APOE4_new',sep = "+"))) #(singlesex delete sex)

######1 COX analysis#####
analyzetable <- filter(analyzetable_backup,AD_day>=1095)
quantile(analyzetable$AD_day)
######1.1 continuous variables#####
var_ofinterest <- temp %>% filter(application_type=='phewas.selected_exposure',calcuType=='linear'|calcuType=='ordered') %>% select(fieldid)
var_ofinterest <- var_ofinterest$fieldid

coeff.old <- data.frame()
for (i in var_ofinterest) {
  analyzetable_noNA <- filter(analyzetable, !is.na(Age_at_baseline),!is.na(sex),!is.na(APOE4_new),
                              !is.na(AD_status),!is.na(AD_day),!is.na(analyzetable[,which(names(analyzetable)==i)])) ##注意读取相应包和文件
  if (nrow(filter(analyzetable_noNA,sex=='male'))==0|nrow(filter(analyzetable_noNA,sex=='female'))==0) {  #for varibale like field 6153 only female 
    FML <- as.formula(paste0('Surv(AD_day,AD_status)~',paste(i,'Age_at_baseline','APOE4_new',sep = "+")))
    pheno_cox_AD <- coxph(FML,data=analyzetable_noNA)
    test.ph <- cox.zph(pheno_cox_AD)
    varid <- i
    coeff <- summary(pheno_cox_AD)$conf.int[1,1]
    lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
    upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
    pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
    singlesex <- 1
    samplesize <- nrow(analyzetable_noNA)
    test.var <- test.ph$table[1,3]
    test.age <- test.ph$table[2,3]
    test.sex <- NA###
    test.apoe4 <- test.ph$table[3,3]
    test.global <- test.ph$table[4,3]
    coeff.new <- data.frame(varid=varid,
                            samplesize=samplesize,
                            coeff=coeff,
                            lower95=lower95,
                            upper95=upper95,
                            pvalue=pvalue,
                            singlesex=singlesex,
                            test.var =test.var,
                            test.age =test.age,
                            test.sex =test.sex,
                            test.apoe4 =test.apoe4,
                            test.global =test.global
    )
    coeff.old <- rbind.data.frame(coeff.old,coeff.new)
    
  } 
  else{ #varibales with both sexes
    FML <- as.formula(paste0('Surv(AD_day,AD_status)~',paste(i,'Age_at_baseline','sex','APOE4_new',sep = "+")))
    pheno_cox_AD <- coxph(FML,data=analyzetable_noNA)
    test.ph <- cox.zph(pheno_cox_AD)
    varid <- i
    coeff <- summary(pheno_cox_AD)$conf.int[1,1]
    lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
    upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
    pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
    singlesex <- 0
    samplesize <- nrow(analyzetable_noNA)
    test.var <- test.ph$table[1,3]
    test.age <- test.ph$table[2,3]
    test.sex <- test.ph$table[3,3]
    test.apoe4 <- test.ph$table[4,3]
    test.global <- test.ph$table[5,3]
    coeff.new <- data.frame(varid=varid,
                            samplesize=samplesize,
                            coeff=coeff,
                            lower95=lower95,
                            upper95=upper95,
                            pvalue=pvalue,
                            singlesex=singlesex,
                            test.var =test.var,
                            test.age =test.age,
                            test.sex =test.sex,
                            test.apoe4 =test.apoe4,
                            test.global =test.global
    )
    
    coeff.old <- rbind.data.frame(coeff.old,coeff.new)
    
  }
}

coeff.1 <- coeff.old
coeff.1

######1.2 binary variables#####
var_ofinterest <- temp %>% filter(application_type=='phewas.selected_exposure',calcuType=='binary') %>% select(fieldid)
var_ofinterest <- var_ofinterest$fieldid

coeff.old <- data.frame()
for (i in var_ofinterest) {
  analyzetable_noNA <- filter(analyzetable, !is.na(Age_at_baseline),!is.na(sex),!is.na(APOE4_new),
                              !is.na(AD_status),!is.na(AD_day),!is.na(analyzetable[,which(names(analyzetable)==i)])) 
  if (nrow(filter(analyzetable_noNA,sex=='male'))==0|nrow(filter(analyzetable_noNA,sex=='female'))==0) {  
    FML <- as.formula(paste0('Surv(AD_day,AD_status)~',paste(i,'Age_at_baseline','APOE4_new',sep = "+")))
    pheno_cox_AD <- coxph(FML,data=analyzetable_noNA)
    test.ph <- cox.zph(pheno_cox_AD)
    varid <- i
    coeff <- summary(pheno_cox_AD)$conf.int[1,1]
    lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
    upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
    pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
    singlesex <- 1
    samplesize <- nrow(analyzetable_noNA)
    test.var <- test.ph$table[1,3]
    test.age <- test.ph$table[2,3]
    test.sex <- NA
    test.apoe4 <- test.ph$table[3,3]
    test.global <- test.ph$table[4,3]
    case_r <- filter(analyzetable_noNA,analyzetable_noNA[,which(names(analyzetable_noNA)==i)]==1)
    control_r <- filter(analyzetable_noNA,analyzetable_noNA[,which(names(analyzetable_noNA)==i)]==0)
    coeff.new <- data.frame(varid=varid,
                            samplesize=paste0(nrow(analyzetable_noNA),'(',
                                              nrow(case_r),'/',
                                              nrow(control_r),')'),
                            coeff=coeff,
                            lower95=lower95,
                            upper95=upper95,
                            pvalue=pvalue,
                            singlesex=singlesex,
                            test.var =test.var,
                            test.age =test.age,
                            test.sex =test.sex,
                            test.apoe4 =test.apoe4,
                            test.global =test.global
    )
    coeff.old <- rbind.data.frame(coeff.old,coeff.new)
    
  } 
  else{ 
    FML <- as.formula(paste0('Surv(AD_day,AD_status)~',paste(i,'Age_at_baseline','sex','APOE4_new',sep = "+")))
    pheno_cox_AD <- coxph(FML,data=analyzetable_noNA)
    test.ph <- cox.zph(pheno_cox_AD)
    varid <- i
    coeff <- summary(pheno_cox_AD)$conf.int[1,1]
    lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
    upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
    pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
    singlesex <- 0
    samplesize <- nrow(analyzetable_noNA)
    test.var <- test.ph$table[1,3]
    test.age <- test.ph$table[2,3]
    test.sex <- test.ph$table[3,3]
    test.apoe4 <- test.ph$table[4,3]
    test.global <- test.ph$table[5,3]
    case_r <- filter(analyzetable_noNA,analyzetable_noNA[,which(names(analyzetable_noNA)==i)]==1)
    control_r <- filter(analyzetable_noNA,analyzetable_noNA[,which(names(analyzetable_noNA)==i)]==0)
    coeff.new <- data.frame(varid=varid,
                            samplesize=paste0(nrow(analyzetable_noNA),'(',
                                              nrow(case_r),'/',
                                              nrow(control_r),')'),
                            coeff=coeff,
                            lower95=lower95,
                            upper95=upper95,
                            pvalue=pvalue,
                            singlesex=singlesex,
                            test.var =test.var,
                            test.age =test.age,
                            test.sex =test.sex,
                            test.apoe4 =test.apoe4,
                            test.global =test.global
    )
    
    coeff.old <- rbind.data.frame(coeff.old,coeff.new)
    
  }
}

coeff.2 <- coeff.old
coeff.2


######1.3 multilevel variables#####
var_ofinterest <- temp %>% filter(application_type=='phewas.selected_exposure',calcuType=='multilevel') %>% select(fieldid)
var_ofinterest <- var_ofinterest$fieldid
var_ofinterest_index <- which(names(analyzetable) %in% var_ofinterest)

coeff.z <- data.frame()
for ( z in var_ofinterest_index) {
  
  analyzetable_noNA <- filter(analyzetable, !is.na(Age_at_baseline),!is.na(sex),!is.na(APOE4_new),
                              !is.na(AD_status),!is.na(AD_day),!is.na(analyzetable[,z])) 
  values0 <- names(table(analyzetable_noNA[,z]))
  coeff.combine <- data.frame()
  for (j in 0:(length(values0)-2)) {
    analyzetable_noNA[,z] <- factor(analyzetable_noNA[,z],levels = values0)
    analyzetable_noNA[,z] <- fct_shift(analyzetable_noNA[,z],j)
    pheno_cox_AD <- coxph(Surv(AD_day,AD_status)~ analyzetable_noNA[,z] + Age_at_baseline + sex+APOE4_new,data=analyzetable_noNA)
    print(summary(pheno_cox_AD))
    test.ph <- cox.zph(pheno_cox_AD)
    coeff.old <- data.frame()
    for ( l in 1:(length(values0)-1-j)) {
      
      varid <- paste0(colnames(analyzetable_noNA)[z],':',j+as.numeric(values0[1])+l,'vs',j+as.numeric(values0[1]))
      coeff <- summary(pheno_cox_AD)$conf.int[l,1]
      lower95 <- summary(pheno_cox_AD)$conf.int[l,3]
      upper95 <- summary(pheno_cox_AD)$conf.int[l,4]
      pvalue <- summary(pheno_cox_AD)$coefficients[l,5]
      singlesex <- NA
      test.var <- test.ph$table[1,3]
      test.age <- test.ph$table[2,3]
      test.sex <- test.ph$table[3,3]
      test.apoe4 <- test.ph$table[4,3]
      test.global <- test.ph$table[5,3]
      
      samplesize <- paste0(nrow(filter(analyzetable_noNA,analyzetable_noNA[,z]==j+as.numeric(values0[1]))),'/',nrow(filter(analyzetable_noNA,analyzetable_noNA[,z]==j+as.numeric(values0[1])+l)))
      coeff.new <- data.frame(varid=varid,
                              samplesize=samplesize,
                              coeff=coeff,
                              lower95=lower95,
                              upper95=upper95,
                              pvalue=pvalue,
                              singlesex=singlesex,
                              test.var =test.var,
                              test.age =test.age,
                              test.sex =test.sex,
                              test.apoe4 =test.apoe4,
                              test.global =test.global)
      
      coeff.old <- rbind.data.frame(coeff.old,coeff.new)
    }
    coeff.combine <- rbind.data.frame(coeff.combine,coeff.old)
  } 
  coeff.z <- rbind.data.frame(coeff.z,coeff.combine)
}
coeff.z
coeff.3 <- coeff.z
coeff.3


######1.4 combine results#####
coeff.f <- rbind.data.frame(coeff.1,coeff.2,coeff.3)
coeff.f.basicCOX.limitto3yearFU <- coeff.f

write.csv(coeff.f.basicCOX.limitto3yearFU,'/AD_phewas/coeff.f.basicCOX.limitto3yearFU.csv',row.names = FALSE)


#####2.COX analysis covar adjusted#####
#'education','BMI','Townsand','Smoking','Alcohal'
analyzetable <- filter(analyzetable_backup,AD_day>=1095) 
analyzetable$Smoking <- factor(analyzetable$Smoking)
analyzetable$Alcohal <- factor(analyzetable$Alcohal)

######2.1.continuous variables#####
var_ofinterest <- temp %>% filter(application_type=='phewas.selected_exposure',calcuType=='linear'|calcuType=='ordered') %>% select(fieldid)
var_ofinterest <- var_ofinterest$fieldid

coeff.old <- data.frame()
for (i in var_ofinterest) {
  analyzetable_noNA <- filter(analyzetable, !is.na(Age_at_baseline),!is.na(sex),!is.na(APOE4_new),
                              !is.na(AD_status),!is.na(AD_day),!is.na(analyzetable[,which(names(analyzetable)==i)]),
                              !is.na(edu_cat),!is.na(Townsand),
                              !is.na(Smoking),!is.na(Alcohal),!is.na(BMI)) 
  if (nrow(filter(analyzetable_noNA,sex=='male'))==0|nrow(filter(analyzetable_noNA,sex=='female'))==0) {  
    FML <- as.formula(paste0('Surv(AD_day,AD_status)~',paste(i,'Age_at_baseline','APOE4_new','edu_cat','Townsand','Smoking','Alcohal','BMI',sep = "+")))
    pheno_cox_AD <- coxph(FML,data=analyzetable_noNA)
    
    varid <- i
    coeff <- summary(pheno_cox_AD)$conf.int[1,1]
    lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
    upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
    pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
    singlesex <- 1
    
    samplesize <- nrow(analyzetable_noNA)
    coeff.new <- data.frame(varid=varid,
                            samplesize=samplesize,
                            coeff=coeff,
                            lower95=lower95,
                            upper95=upper95,
                            pvalue=pvalue,
                            singlesex=singlesex
    )
    coeff.old <- rbind.data.frame(coeff.old,coeff.new)
    
  } 
  else{ 
    FML <- as.formula(paste0('Surv(AD_day,AD_status)~',paste(i,'Age_at_baseline','sex','APOE4_new','edu_cat','Townsand','Smoking','Alcohal','BMI',sep = "+")))
    pheno_cox_AD <- coxph(FML,data=analyzetable_noNA)
    
    varid <- i
    coeff <- summary(pheno_cox_AD)$conf.int[1,1]
    lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
    upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
    pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
    singlesex <- 0
    samplesize <- nrow(analyzetable_noNA)
    coeff.new <- data.frame(varid=varid,
                            samplesize=samplesize,
                            coeff=coeff,
                            lower95=lower95,
                            upper95=upper95,
                            pvalue=pvalue,
                            singlesex=singlesex
    )
    
    coeff.old <- rbind.data.frame(coeff.old,coeff.new)
    
  }
}
coeff.1 <- coeff.old
coeff.1

######2.2.binary variables#####
var_ofinterest <- temp %>% filter(application_type=='phewas.selected_exposure',calcuType=='binary') %>% select(fieldid)
var_ofinterest <- var_ofinterest$fieldid

coeff.old <- data.frame()
for (i in var_ofinterest) {
  analyzetable_noNA <- filter(analyzetable, !is.na(Age_at_baseline),!is.na(sex),!is.na(APOE4_new),
                              !is.na(AD_status),!is.na(AD_day),!is.na(analyzetable[,which(names(analyzetable)==i)]),
                              !is.na(edu_cat),is.na(Townsand),
                              !is.na(Smoking),!is.na(Alcohal),!is.na(BMI)) 
  if (nrow(filter(analyzetable_noNA,sex=='male'))==0|nrow(filter(analyzetable_noNA,sex=='female'))==0) {  
    FML <- as.formula(paste0('Surv(AD_day,AD_status)~',paste(i,'Age_at_baseline','APOE4_new','edu_cat','Townsand','Smoking','Alcohal','BMI',sep = "+")))
    pheno_cox_AD <- coxph(FML,data=analyzetable_noNA)
    
    varid <- i
    coeff <- summary(pheno_cox_AD)$conf.int[1,1]
    lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
    upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
    pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
    singlesex <- 1
    samplesize <- nrow(analyzetable_noNA)
    
    case_r <- filter(analyzetable_noNA,analyzetable_noNA[,which(names(analyzetable_noNA)==i)]==1)
    control_r <- filter(analyzetable_noNA,analyzetable_noNA[,which(names(analyzetable_noNA)==i)]==0)
    coeff.new <- data.frame(varid=varid,
                            samplesize=paste0(nrow(analyzetable_noNA),'(',
                                              nrow(case_r),'/',
                                              nrow(control_r),')'),
                            coeff=coeff,
                            lower95=lower95,
                            upper95=upper95,
                            pvalue=pvalue,
                            singlesex=singlesex
    )
    coeff.old <- rbind.data.frame(coeff.old,coeff.new)
    
  } 
  else{ #两个性别都有变量分布
    FML <- as.formula(paste0('Surv(AD_day,AD_status)~',paste(i,'Age_at_baseline','sex','APOE4_new','edu_cat','Townsand','Smoking','Alcohal','BMI',sep = "+")))
    pheno_cox_AD <- coxph(FML,data=analyzetable_noNA)
    varid <- i
    
    coeff <- summary(pheno_cox_AD)$conf.int[1,1]
    lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
    upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
    pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
    singlesex <- 0
    samplesize <- nrow(analyzetable_noNA)
    
    case_r <- filter(analyzetable_noNA,analyzetable_noNA[,which(names(analyzetable_noNA)==i)]==1)
    control_r <- filter(analyzetable_noNA,analyzetable_noNA[,which(names(analyzetable_noNA)==i)]==0)
    coeff.new <- data.frame(varid=varid,
                            samplesize=paste0(nrow(analyzetable_noNA),'(',
                                              nrow(case_r),'/',
                                              nrow(control_r),')'),
                            coeff=coeff,
                            lower95=lower95,
                            upper95=upper95,
                            pvalue=pvalue,
                            singlesex=singlesex
                            
    )
    
    coeff.old <- rbind.data.frame(coeff.old,coeff.new)
    
  }
}

coeff.2 <- coeff.old
coeff.2


######2.3.multilevel variables#####
var_ofinterest <- temp %>% filter(application_type=='phewas.selected_exposure',calcuType=='multilevel') %>% select(fieldid)
var_ofinterest <- var_ofinterest$fieldid
var_ofinterest_index <- which(names(analyzetable) %in% var_ofinterest)

coeff.z <- data.frame()
for ( z in var_ofinterest_index) {
  
  analyzetable_noNA <- filter(analyzetable, !is.na(Age_at_baseline),!is.na(sex),!is.na(APOE4_new),
                              !is.na(AD_status),!is.na(AD_day),!is.na(analyzetable[,z]),
                              #!is.na(X6138.0.0)&X6138.0.0!=-3,
                              !is.na(edu_cat),!is.na(Townsand),
                              !is.na(Smoking),!is.na(Alcohal),!is.na(BMI)) ##注意读取相应包和文件
  values0 <- names(table(analyzetable_noNA[,z]))###也可以这样： values0 <- unique(analyzetable_noNA[,z])
  coeff.combine <- data.frame()
  for (j in 0:(length(values0)-2)) {
    analyzetable_noNA[,z] <- factor(analyzetable_noNA[,z],levels = values0)
    analyzetable_noNA[,z] <- fct_shift(analyzetable_noNA[,z],j)
    pheno_cox_AD <- coxph(Surv(AD_day,AD_status)~ analyzetable_noNA[,z] + Age_at_baseline + sex+APOE4_new+edu_cat+Townsand+Smoking+Alcohal+BMI,data=analyzetable_noNA)
    
    coeff.old <- data.frame()
    for ( l in 1:(length(values0)-1-j)) {
      
      varid <- paste0(colnames(analyzetable_noNA)[z],':',j+as.numeric(values0[1])+l,'vs',j+as.numeric(values0[1]))
      coeff <- summary(pheno_cox_AD)$conf.int[l,1]
      lower95 <- summary(pheno_cox_AD)$conf.int[l,3]
      upper95 <- summary(pheno_cox_AD)$conf.int[l,4]
      pvalue <- summary(pheno_cox_AD)$coefficients[l,5]
      singlesex <- NA
      
      
      samplesize <- paste0(nrow(filter(analyzetable_noNA,analyzetable_noNA[,z]==j+as.numeric(values0[1])+l)),'/',nrow(filter(analyzetable_noNA,analyzetable_noNA[,z]==j+as.numeric(values0[1]))))
      coeff.new <- data.frame(varid=varid,
                              samplesize=samplesize,
                              coeff=coeff,
                              lower95=lower95,
                              upper95=upper95,
                              pvalue=pvalue,
                              singlesex=singlesex
      )
      
      coeff.old <- rbind.data.frame(coeff.old,coeff.new)
    }
    coeff.combine <- rbind.data.frame(coeff.combine,coeff.old)
  } 
  coeff.z <- rbind.data.frame(coeff.z,coeff.combine)
}
coeff.z
coeff.3 <- coeff.z
coeff.3

######2.4.combine results#####
coeff.f <- rbind.data.frame(coeff.1,coeff.2,coeff.3)
coeff.f.crct <- coeff.f

write.csv(coeff.f.educrct,'/AD_phewas/coeff.f.crct.csv',row.names = FALSE)

#============================================================================#
#============================================================================#




