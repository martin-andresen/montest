setwd("C:/Users/martiea/Dropbox/Prosjekter/working/testing IV monotonicity/Tymon's non-Young replications/completed replications")
library(haven)
library(data.table)
library(grf)
library(fixest)
library(tinyplot)

montests=list()

#ALLCOTT 2020
data=data.table(read_dta("Allcott2020_data.dta"))
data=data[,c("v_wins","stratum","T","D","endline_wta_update")]
montests[["Allcott2020"]]=montest(data=data,Y="endline_wta_update",X=c("v_wins","stratum"),Z="T",D="D",test="simple")
#ivreg2 endline_wta_update v_wins i.stratum (D = T), robust
#NOTE: ONLY ONE REGRESSION IN THIS REPLICATION FILE?
#MAIN SPEC IS FROM APPENDIX TABLE A12

#AMBRUS 2020: MULTIPLE INSTRUMENTS
#ivreg2 log_rentals_1864 (death_ind = broad##c.dist_broad) dist_urinal old_sewer dist_pump if dist_netw<=0.292, cl(block)

#ASHER 2020
##MAIN EST: TABLE 3
##ivregress 2sls transport_index_andrsn (r2012 = t) left right primary_school med_center elect tdist irr_share ln_land pc01_lit_share pc01_sc_share bpl_landed_share bpl_inc_source_sub_share bpl_inc_250plus i.vhg_dist_id [aw = kernel_tri_ik] if mainsample, vce(robust)
data=data.table(read_dta("Asher2020_data.dta"))
data=data[mainsample==1]
data[,runvar:=ifelse(t==1,right,left)]
data=data[,c("transport_index_andrsn","r2012","t","runvar","primary_school","med_center","elect","tdist","irr_share","ln_land","pc01_lit_share","pc01_sc_share","bpl_landed_share","bpl_inc_source_sub_share","bpl_inc_250plus","vhg_dist_id","kernel_tri_ik")]
montests[["Asher2020"]]=montest(data=data,Z="t",D="r2012",X=c("runvar","primary_school","med_center","elect","tdist","irr_share","ln_land","pc01_lit_share","pc01_sc_share","bpl_landed_share","bpl_inc_source_sub_share","bpl_inc_250plus","vhg_dist_id"),test="simple")
##NOT NP ID - fuzzy RD

##AUTOR 2020a
#ivreg2 dhs2_tot_cont_2002_2010 (d_imp_usch_pd=d_imp_otch_lag_pd) [aw=sh_district_2002], cluster(czone congressionaldistrict)
data=data.table(read_dta("Autor2020a_data_1.dta"))
X=c("reg_encen","reg_midatl","reg_wncen","reg_satl","reg_escen","reg_wscen","reg_mount","reg_pacif","l_shind_manuf_cbp","l_sh_routine33","l_task_outsource","shnr_pres2000","shnr_pres1996","l_sh_pop_f","l_sh_pop_edu_c","l_sh_fborn","l_sh_pop_age_1019","l_sh_pop_age_2029","l_sh_pop_age_3039","l_sh_pop_age_4049","l_sh_pop_age_5059","l_sh_pop_age_6069","l_sh_pop_age_7079","l_sh_pop_age_8000","l_sh_pop_white","l_sh_pop_black","l_sh_pop_asian","l_sh_pop_hispanic")
cols=c(X,"d_imp_usch_pd","d_imp_otch_lag_pd","czone","sh_district_2002")
data=data[, ..cols]
montests[["Autor2020a"]]=montest(data=data,Z="d_imp_otch_lag_pd",D="d_imp_usch_pd",X=X,test="simple",weight="sh_district_2002",cluster="czone")
#MAIN SPEC TABLE 3, row 1, col 5
#ivreg2 dhs2_tot_cont_2002_2010 (d_imp_usch_pd=d_imp_otch_lag_pd) reg* l_shind_manuf_cbp l_sh_routine33 l_task_outsource shnr_pres2000 shnr_pres1996 l_sh_pop_f l_sh_pop_edu_c l_sh_fborn l_sh_pop_age_1019 l_sh_pop_age_2029 l_sh_pop_age_3039 l_sh_pop_age_4049 l_sh_pop_age_5059 l_sh_pop_age_6069 l_sh_pop_age_7079 l_sh_pop_age_8000 l_sh_pop_white l_sh_pop_black l_sh_pop_asian l_sh_pop_hispanic [aw=sh_district_2002], cluster(czone congressionaldistrict)
##OOPS: TWO-WAY CLUSTERING, imp??lemented cluster on czone only!


##BANDIERA 2020
#ivregress 2sls Qcontrol_body control_body _B* age (QC_clubparticipateIMP=treatment) if panel==1, cluster(villid) first
data=data.table(read_dta("Bandiera2020_data.dta"))
data[,Bbranch := as.integer(sub("^_Bbranch_na_", "", branch_cols[max.col(.SD)])),.SDcols = branch_cols]
X=c("control_body reg_encen","reg_midatl","reg_wncen","reg_satl","reg_escen","reg_wscen","reg_mount","reg_pacif","l_shind_manuf_cbp","l_sh_routine33","l_task_outsource","shnr_pres2000","shnr_pres1996","l_sh_pop_f","l_sh_pop_edu_c","l_sh_fborn","l_sh_pop_age_1019","l_sh_pop_age_2029","l_sh_pop_age_3039","l_sh_pop_age_4049","l_sh_pop_age_5059","l_sh_pop_age_6069","l_sh_pop_age_7079","l_sh_pop_age_8000","l_sh_pop_white","l_sh_pop_black","l_sh_pop_asian","l_sh_pop_hispanic")
cols=c(X,"d_imp_usch_pd","d_imp_otch_lag_pd","czone","sh_district_2002")
data=data[, ..cols]

##BURSZTYN 2020
data=data.table(read_dta("Bursztyn2020_data_1.dta"))
data=data[data$public==1]
data[,sob_culture_50 := sob_culture - 50]
data=data[,c("female","age","married","years_edu","income_000s","white","donate","trump","sob_culture_50")]
montests[["Bursztyn2020"]]=montest(data=data,Z="trump",D="sob_culture_50",X=c("female","age","married","years_edu","income_000s","white"),test="simple")
#ivregress 2sls donate (sob_culture_50 = trump) female age married years_edu income_000s white
#MAIN estimate: table 2, panel B, col 6.



##CAPRETTINI 2020
data=data.table(read_dta("Caprettini2020_data.dta"))
data=data[sample==1]
data=data[,c("thresh","heavysh","cer","log_density","agri_share","log_sex_ratio","log_distel","log_distnews","REGION")]
montests[["Caprettini2020"]]=montest(data=data,D="thresh",Z="heavysh",X=c("cer","log_density","agri_share","log_sex_ratio","log_distel","log_distnews","REGION"),test="simple")
#ivreg2 SWING (thresh = heavysh) cer log_density agri_share log_sex_ratio log_distel log_distnews _IREGION_* if sample == 1, r first
#MAIN ESTIMATE: table 2, col 8

data=data.table(read_dta("Dupas2013_data.dta"))
data=data[,c("bank_savings","active","treatment","wave2","wave3","bg_boda","bg_malevendor","bg_boda_wave2","bg_malevendor_wave2","bg_married","bg_num_children","bg_age","bg_kis_read","bg_rosca_contrib_lyr","filled_log")]
montests[["Dupas2013"]]=montest(data=data,X=c("wave2","wave3","bg_boda","bg_malevendor","bg_boda_wave2","bg_malevendor_wave2","bg_married","bg_num_children","bg_age","bg_kis_read","bg_rosca_contrib_lyr","filled_log"),Y="bank_savings",D="active",Z="treatment",test="simple")
#ivreg bank_savings (active = treatment) wave2 wave3 bg_boda bg_malevendor bg_boda_wave2 bg_malevendor_wave2 bg_married bg_num_children bg_age bg_kis_read bg_rosca_contrib_lyr filled_log, robust

data=data.table(read_dta("Galiani2011_data_1.dta"))
data=data[cohort>1957&cohort<1963]
data=data[,c("crimerate","sm","highnumber","cohort")]
montests[["Galiani2011"]]=montest(data=data,X="cohort",Y="crimerate",D="sm",Z="highnumber",test="simple")
#i: ivreg crimerate (sm = highnumber) i.cohort if cohort > 1957 & cohort < 1963, robust

data=data.table(read_dta("Gerber2020_data.dta")) # cLOSE TO
data=data[,c("votemarg_post","t_close","ppstaten","vote_admin2000","vote_admin2002","vote_admin2004","vote_admin2006","vote_admin2008")]
montests[["Gerber2020"]]=montest(data=data,,D="votemarg_post",Z="t_close",X=c("ppstaten","vote_admin2000","vote_admin2002","vote_admin2004","vote_admin2006","vote_admin2008"),test="simple")
#xi: ivreg2 vote_admin1 (votemarg_post=t_close) i.ppstaten vote_admin2000 vote_admin2002 vote_admin2004 vote_admin2006 vote_admin2008, first savefirst robust

data=data.table(read_dta("Gregg2020_data.dta")) ## NO WORK
data=data[,c("Form","RelRevCPMinus","Industry","PIHerfIndex")]
data[,Industry:=as.numeric(as.factor(Industry))]
data=na.omit(data)
#montests[["Gregg2020"]]=montest(data=data,D="Form",Z="RelRevCPMinus",X=c("Industry","PIHerfIndex"),test="simple")
#xi: ivregress 2sls logRevperWorker (Form = RelRevCPMinus) i.Industry PIHerfIndex, r first

data=data.table(read_dta("Guryan2010_data_1.dta")) ##REJECT, but because of minsize bug
data[, dyrmo := {
  M   <- as.matrix(.SD) * 1L              # coerce logicals to 0/1 if needed
  rs  <- rowSums(M)                        # count of 1s per row
  idx <- max.col(M, ties.method = "first") # position of the 1
  idx[rs != 1L] <- NA_integer_             # require exactly one 1
  labs <- sub("^dyrmo[_]?", "", names(.SD))# labels from colnames
  factor(labs[idx], levels = labs)
}, .SDcols = patterns("^dyrmo")]
data=data[,c("lzsales1","gzanywin","lzsales","dyrmo")]
data=na.omit(data)
data[,dyrmo:=as.numeric(dyrmo)]
montests[["Guryan2010"]]=montest(data=data,D="lzsales1",Z="gzanywin",X=c("lzsales","dyrmo"),test="simple")
#ivregress 2sls lzsales2 (lzsales1 = gzanywin) lzsales dyrmo*, robust first

#data=data.table(read_dta("Lagos2020_data.dta")) #NO X
#data=data[(FOMC_Hbased==1|FOMC_Hbased==0)&date>="1994-01-01"&date<="2008-01-01"]
#data=data[,c("dr","wr")]
#ivregress 2sls ret_m (dr=wr) if (FOMC_Hbased==1|FOMC_Hbased==0)&date>=d(01jan1994)&date<=d(01jan2008), robust

data=data.table(read_dta("Owens2020_data_1.dta"))
data=data[,c("log_dist_to_park","log_dist_to_highway","log_dist_to_airport","log_dist_to_water","log_dist_to_college","logRj","logAi")]
montests[["Owens2020"]]=montest(data=data,X=c("log_dist_to_park","log_dist_to_highway","log_dist_to_airport","log_dist_to_water","log_dist_to_college"),D="logRj",Z="logAi",test="simple")
#ivreg2 logamenities log_dist_to_park log_dist_to_highway log_dist_to_airport log_dist_to_water log_dist_to_college (logRj = logAi), robust

data=data.table(read_dta("Pons2018_data.dta"))
data=data[merge_results12==1]
data=data[,c("allocated","treatment","stratum_identifier")]
montests[["Pons2018"]]=montest(data=data,X="stratum_identifier",D="allocated",Z="treatment",test="simple")
#xtivreg2 prop_turnout_pr12t1_an (allocated = treatment) if merge_results12 == 1, i(stratum_identifier) fe robust
