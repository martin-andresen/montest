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
feols(D~T+v_wins +i(stratum),data=data) ##FS
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
feols(r2012~t+left+right+primary_school+med_center+elect+tdist+irr_share+ln_land+pc01_lit_share+bpl_landed_share+bpl_inc_source_sub_share+bpl_inc_250plus+i(vhg_dist_id),data=data,weight=~kernel_tri_ik)
X=c("runvar","primary_school","med_center","elect","tdist","irr_share","ln_land","pc01_lit_share","pc01_sc_share","bpl_landed_share","bpl_inc_source_sub_share","bpl_inc_250plus","vhg_dist_id")
cols=c(X,"transport_index_andrsn","r2012","t","kernel_tri_ik")
data=data[,..cols]
##montests[["Asher2020"]]=montest(data=data,Z="t",D="r2012",X=X,test="simple")
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
##NOT NP ID?

##BANDIERA 2020
#ivregress 2sls Qcontrol_body control_body _B* age (QC_clubparticipateIMP=treatment) if panel==1, cluster(villid) first
data=data.table(read_dta("Bandiera2020_data.dta"))
data=data[panel==1]
X=c("control_body","branchno","age")
cols=c(X,"treatment","QC_clubparticipateIMP","villid")
feols(QC_clubparticipateIMP~treatment+control_body+i(branchno)+age, data=data)
data=data[,..cols]
montests[["Bandiera2020"]]=montest(data=data,Z="treatment",D="QC_clubparticipateIMP",X=X,test="simple",cluster="villid")


##Bau2020
#ivreg2 TVA_mean (mean_1 = mean_2) female local some_training BA_plus lessthan4 temp_contract _Idistrict__* if gov==1, cluster(group) endog(mean_1)
data=data.table(read_dta("Bau2020_data_1.dta"))
data[,district:=as.integer(as.factor(district_name))]
data=data[gov==1]
X=c("schoolid","female","local","some_training","BA_plus","lessthan4","temp_contract")
cols=c(X,"mean_1","mean_2","group","district")
data=data[,..cols]
feols(as.formula(paste0(c("mean_1~mean_2",X,"i(district)"),collapse="+")),data=data)
X=c(X,"district")
montests[["Bau2020"]]=montest(data=data,D="mean_1",Z="mean_2",X=X,test="simple",cluster="group")
##MAIN ESTIMATE: table 3, cols 6-7 (district or school id)

##Becker2019
#ivreg2 share_antisem_reich_1890 (f_prot_1882=kmwittenberg) f_young f_fem f_ortsgeb f_pruss hhsize lnpop posen f_urban, cluster(code_reichstag_wk)
#Main spec: table 4, panel C. However, table 6 might be thought of as the main est, we don't ahve data.
data=data.table(read_dta("Becker2019_data.dta"))
X=c("f_young","f_fem","f_ortsgeb","f_pruss","hhsize","lnpop","posen","f_urban")
cols=c(X,"f_prot_1882","kmwittenberg","code_reichstag_wk")
data=data[,..cols]
feols(as.formula(paste0(c("f_prot_1882~kmwittenberg",X),collapse="+")),data=data) ##NEGATIVE FS
data[,kmwittenberg:=-kmwittenberg]
montests[["Becker2019"]]=montest(data=data,D="f_prot_1882",Z="kmwittenberg",X=X,test="simple",cluster="code_reichstag_wk")

#Bergqquist2020
#ivregress 2sls weighted_price_adj_trim _Iweek_* _Imarket_na_* (num_entrants = S2) if S1!=1 [aweight=num_traders_inv], r cluster(market_block)
data=data.table(read_dta("Bergquist2020_data.dta"))
data[,marketno:=as.factor(market_name)]
data[,marketno:=as.numeric(marketno)]
data[,market_block:=as.numeric(as.factor(market_block))]
X=c("week","marketno")
data=data[S1!=1]
cols=c(X,"weighted_price_adj_trim","num_entrants","S2","num_traders_inv","market_block")
data=data[,..cols]
feols(num_entrants~S2+i(marketno)+i(week),weight=~num_traders_inv,data=data)
montests[["Bergquist2020"]]=montest(data=data,D="num_entrants",Z="S2",X=X,test="simple",weight="num_traders_inv",cluster="market_block")
#Main estimate: table 5, col 4

##Bound2020
#xtivreg2 l_foreign_fresh (l_state_ap = l_state_app_at_state) l_population y1-y17 [w=weight] if Research==1, fe cluster(state_of_college)
data=data.table(read_dta("Bound2020_data_1.dta"))
data=data[Research==1]
yvars <- paste0("y", 1:17)
data[, year := max.col(as.matrix(.SD)), .SDcols = yvars]
X=c("l_population","unitid","year")
cols=c(X,"weight","l_state_ap","l_state_app_at_state","l_foreign_fresh","state_of_college")
data=data[,..cols]
data[,state_of_college:=as.numeric(as.factor(state_of_college))]
montests[["Bound2020"]]=montest(data=data,D="l_state_ap",Z="l_state_app_at_state",X=X,test="simple",weight="weight")
##MAIN ESTIMATE: Table 2, col2
##OBS cant cluster by state - too few!



##BURSZTYN 2020
data=data.table(read_dta("Bursztyn2020_data_1.dta"))
data=data[data$public==1]
data[,sob_culture_50 := sob_culture - 50]
X=c("female","age","married","years_edu","income_000s","white")
cols=c(X,"donate","trump","sob_culture_50")
data=data[,..cols]
montests[["Bursztyn2020"]]=montest(data=data,Z="trump",D="sob_culture_50",X=X,test="simple")
feols(as.formula(paste0(c("sob_culture_50~trump",X),collapse="+")),data=data)
#ivregress 2sls donate (sob_culture_50 = trump) female age married years_edu income_000s white
#MAIN estimate: table 2, panel B, col 6.

##Butters2020
#ivregress 2sls logcaputil logwage logrent logelecprice logpipc logpi logempllh ur logllhshare logllhsqmile logpipc5yr logpi5yr logempllh5yr i.segment i.year (demandvolatility = instrument1_imt), vce(cluster metro) first
#main spec: table 2 col c
data=data.table(read_dta("Butters2020_data_1.dta"))
X=c("logwage", "logrent", "logelecprice", "logpipc", "logpi", "logempllh",
    "ur", "logllhshare", "logllhsqmile",
    "logpipc5yr", "logpi5yr", "logempllh5yr")
feols(as.formula(paste0(c("demandvolatility~instrument1_imt",X,"i(segment)","i(year)"),collapse="+")),data=data)
X=c(X,"segment", "year")
cols=c(X,"logcaputil","demandvolatility","instrument1_imt","metro")
data=data[,..cols]
montests[["Butters2020"]]=montest(data=data,Z="instrument1_imt",D="demandvolatility",X=X,test="simple")
#NOTE: variable "intstrument" doesn't exist - assumed instrment1_imt
##couldn't cluster - 92 clusters


##Cahyadi2020
#ivregress 2sls pre_natal_visits hh_head_agr_baseline_nm hh_head_serv_baseline_nm hh_educ*baseline_nm roof_type*baseline_nm wall_type*baseline_nm floor_type*baseline_nm clean_water_baseline_nm own_latrine_baseline_nm square_latrine_baseline_nm own_septic_tank_baseline_nm electricity_PLN_baseline_nm hhsize_ln_baseline_nm logpcexp_baseline_nm *miss kabu_* (pkh_by_this_wave = L07) if survey_round == 2, vce(cluster kecamatan)
#Main spec: table 2, col 2
data=data.table(read_dta("Cahyadi2020_data_1.dta"))
data=data[survey_round==2]

vars <- paste0("kabu_", 1:28)
data[, kabu := {
  m <- as.matrix(.SD)
  out <- max.col(m, ties.method = "first")
  out[is.na(rowSums(m)) == TRUE] <- 0
  out
}, .SDcols = vars]

vars <- paste0("hh_educ_", 1:9,"_baseline_nm")
data[, hh_educ := {
  m <- as.matrix(.SD)
  out <- max.col(m, ties.method = "first")
  out[is.na(rowSums(m)) == TRUE] <- 0
  out
}, .SDcols = vars]

vars <- paste0("floor_type", 1:9,"_baseline_nm")
data[, floor_type := {
  m <- as.matrix(.SD)
  out <- max.col(m, ties.method = "first")
  out[is.na(rowSums(m)) == TRUE] <- 0
  out
}, .SDcols = vars]

vars <- paste0("roof_type", 1:9,"_baseline_nm")
data[, roof_type := {
  m <- as.matrix(.SD)
  out <- max.col(m, ties.method = "first")
  out[is.na(rowSums(m)) == TRUE] <- 0
  out
}, .SDcols = vars]

vars <- paste0("wall_type", 1:9,"_baseline_nm")
data[, wall_type := {
  m <- as.matrix(.SD)
  out <- max.col(m, ties.method = "first")
  out[is.na(rowSums(m)) == TRUE] <- 0
  out
}, .SDcols = vars]

X=c("hh_head_agr_baseline_nm", "hh_head_serv_baseline_nm", "clean_water_baseline_nm", "own_latrine_baseline_nm", "square_latrine_baseline_nm", "own_septic_tank_baseline_nm", "electricity_PLN_baseline_nm", "hhsize_ln_baseline_nm", "logpcexp_baseline_nm"
,"hh_educ","clean_water_baseline_miss","own_latrine_baseline_miss","square_latrine_baseline_miss","own_septic_tank_baseline_miss","electricity_PLN_baseline_miss","logpcexp_baseline_miss","hhsize_ln_baseline_miss")
feols(as.formula(paste0(c("pkh_by_this_wave~L07",X,"i(kabu)","i(floor_type)","i(roof_type)","i(wall_type)"),collapse="+")),data=data)
X=c(X,"kabu","floor_type","wall_type","roof_type")

cols=c(X,"pre_natal_visits","pkh_by_this_wave","L07","kecamatan")
data=data[,..cols]
montests[["Cahyadi2020"]]=montest(data=data,Z="L07",D="pkh_by_this_wave",X=X,test="simple")

#Cai 2015
#ivregress 2sls takeup_survey (pre_takeup_rate = default) male age agpop ricearea_2010 literacy intensive risk_averse disaster_prob _Ivilid* if delay == 1 & info_none == 0, cluster(address)
##baseline spec: table 6, col 5 - oops should probably be col 3, but can't find the spec
data=data.table(read_dta("Cai2015_data_1.dta"))
X=c("male", "age", "agpop", "ricearea_2010", "literacy", "intensive", "risk_averse", "disaster_prob")
feols(as.formula(paste0(c("pre_takeup_rate~default",X,"i(vilid)"),collapse="+")),data=data)
cols=c(X,"pre_takeup_rate","takeup_survey","takeup_survey","address","default")
data=data[,..cols]
data[,address:=as.integer(as.factor(address))]
montests[["Cai2015"]]=montest(data=data,Z="default",D="pre_takeup_rate",X=X,test="simple",cluster="address")

##CAPRETTINI 2020
data=data.table(read_dta("Caprettini2020_data.dta"))
data=data[sample==1]
X=c("cer","log_density","agri_share","log_sex_ratio","log_distel","log_distnews")
cols=c(X,"thresh","heavysh","SWING","REGION")
data=data[,..cols]
feols(as.formula(paste0(c("thresh~heavysh",X,"i(REGION)"),collapse="+")),data=data) ##NEGATIVE FS
data[,heavysh:=-heavysh]
X=c(X,"REGION")
montests[["Caprettini2020"]]=montest(data=data,D="thresh",Z="heavysh",X=X,test="simple")
#ivreg2 SWING (thresh = heavysh) cer log_density agri_share log_sex_ratio log_distel log_distnews _IREGION_* if sample == 1, r first
#MAIN ESTIMATE: table 2, col 81

##Carrol 2020: Multiple instruments

##Collins 2013
#ivreg2 lnfaminc80 pownocc50 lnmedval50 pdilap50 poldunits50 punitswoplumb50 pcrowd50 lnpop50 pnonwht50 plf_manuf_50 pemp50 medsch50 lnfaminc50 pinc_under2g_50 _I* (app_funds_pc50=yrsexposure_UR), robust cluster(statefip)
#baseline: table 3 panel A
data=data.table(read_dta("Collins2013_data_1.dta"))
X=c("pownocc50", "lnmedval50", "pdilap50", "poldunits50", "punitswoplumb50",
    "pcrowd50", "lnpop50", "pnonwht50", "plf_manuf_50", "pemp50",
    "medsch50", "lnfaminc50", "pinc_under2g_50")
cols=c(X,"region","lnfaminc80","app_funds_pc50","yrsexposure_UR","statefip")
data=data[,..cols]
feols(as.formula(paste0(c("app_funds_pc50~yrsexposure_UR",X,"i(region)"),collapse="+")),data=data) ##NEGATIVE FS
montests[["Collins2013"]]=montest(data=data,D="app_funds_pc50",Z="yrsexposure_UR",X=X,test="simple")

##DeMel 2013
#xi: ivreg2 realprofits (registerDS= treat2 treat3 treat4) lagrealprofits i.strata wave3-wave8 if treat1==0, robust cluster(sheno)
#Multiple instruments
#baseline table 5A - col 1

#Draca 2011
#xi: ivreg dltot (dlhpop=full_treat1) full i.week [aw=pop], cluster(ocu)
#Baseline: Table 2 panel C) col 3
data=data.table(read_dta("Draca2011_data_1.dta"))
X=c("full","week")
cols=c(X,"dlhpop","full_treat1","pop","ocu")
data=data[,..cols]
feols(dlhpop~full_treat1+full+i(week),data=data,weight=~pop,cluster=~ocu) ##NEGATIVE FS
montests[["Draca2011"]]=montest(data=data,D="dlhpop",Z="full_treat1",X=X,test="simple",weight="pop")
##Couldn't cluster!

##Dupas 2013
data=data.table(read_dta("Dupas2013_data.dta"))
data=data[,c("bank_savings","active","treatment","wave2","wave3","bg_boda","bg_malevendor","bg_boda_wave2","bg_malevendor_wave2","bg_married","bg_num_children","bg_age","bg_kis_read","bg_rosca_contrib_lyr","filled_log")]
X=c("wave2","wave3","bg_boda","bg_malevendor","bg_boda_wave2","bg_malevendor_wave2","bg_married","bg_num_children","bg_age","bg_kis_read","bg_rosca_contrib_lyr","filled_log")
cols=c(X,"bank_savings","active","treatment")
data=data[,..cols]
feols(as.formula(paste0(c("active~treatment",X),collapse="+")),data=data)
montests[["Dupas2013"]]=montest(data=data,X=X,D="active",Z="treatment",test="simple")
#ivreg bank_savings (active = treatment) wave2 wave3 bg_boda bg_malevendor bg_boda_wave2 bg_malevendor_wave2 bg_married bg_num_children bg_age bg_kis_read bg_rosca_contrib_lyr filled_log, robust

##Galiani 2011
data=data.table(read_dta("Galiani2011_data_1.dta"))
data=data[cohort>1957&cohort<1963]
data=data[,c("crimerate","sm","highnumber","cohort")]
feols(sm~highnumber+i(cohort),data=data)
montests[["Galiani2011"]]=montest(data=data,X="cohort",D="sm",Z="highnumber",test="simple")
#i: ivreg crimerate (sm = highnumber) i.cohort if cohort > 1957 & cohort < 1963, robust


##Gerber 2020
data=data.table(read_dta("Gerber2020_data.dta"))
data=data[,c("votemarg_post","t_close","ppstaten","vote_admin2000","vote_admin2002","vote_admin2004","vote_admin2006","vote_admin2008")]
X=c("ppstaten","vote_admin2000","vote_admin2002","vote_admin2004","vote_admin2006","vote_admin2008")
feols(votemarg_post~t_close+i(ppstaten)+vote_admin2002+vote_admin2004+vote_admin2008,data=data)
data[,t_close:=-t_close] ##neagtive FS
montests[["Gerber2020"]]=montest(data=data,D="votemarg_post",Z="t_close",X=X,test="simple")
#xi: ivreg2 vote_admin1 (votemarg_post=t_close) i.ppstaten vote_admin2000 vote_admin2002 vote_admin2004 vote_admin2006 vote_admin2008, first savefirst robust
##Not np identified - fuzzy RD?



#Glitz 2020
#ivreg2 c3difflnTFP (espionage=exp_inf_gva_old2) difflnTFP patents br_* yd_* [aw=weight_workers], cluster(branch)
data=data.table(read_dta("Glitz2020_data_1.dta"))
vars <- paste0("yd_", 1:18)
data[, yd := {
  m <- as.matrix(.SD)
  out <- max.col(m, ties.method = "first")
  out[is.na(rowSums(m)) == TRUE] <- 0
  out
}, .SDcols = vars]
vars <- paste0("br_", 1:16)
data[,  br := {
  m <- as.matrix(.SD)
  out <- max.col(m, ties.method = "first")
  out[is.na(rowSums(m)) == TRUE] <- 0
  out
}, .SDcols = vars]
X=c("difflnTFP","patents","br","yd")
cols=c(X,"espionage","exp_inf_gva_old2","weight_workers","branch","c3difflnTFP")
data=data[,..cols]
feols(espionage~exp_inf_gva_old2+patents+difflnTFP+i(yd)+i(br),weight=~weight_workers,cluster=~branch,data=data)
montests[["Glitz2020"]]=montest(data=data,D="espionage",Z="exp_inf_gva_old2",X=X,weight="weight_workers",test="simple")
##COULDN'T cluster

#Gregg2020
data=data.table(read_dta("Gregg2020_data.dta")) ## NO WORK
data=data[,c("Form","RelRevCPMinus","Industry","PIHerfIndex")]
data[,Industry:=as.numeric(as.factor(Industry))]
data=na.omit(data)
feols(Form~RelRevCPMinus+i(Industry)+PIHerfIndex,data=data)
montests[["Gregg2020"]]=montest(data=data,D="Form",Z="RelRevCPMinus",X=c("Industry","PIHerfIndex"),test="simple")
#xi: ivregress 2sls logRevperWorker (Form = RelRevCPMinus) i.Industry PIHerfIndex, r first

##Guryan 2010
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
data[,dyrmo:=as.integer(dyrmo)]
feols(lzsales1~gzanywin+lzsales+i(dyrmo),data=data)
montests[["Guryan2010"]]=montest(data=data,D="lzsales1",Z="gzanywin",X=c("lzsales","dyrmo"),test="simple")
#ivregress 2sls lzsales2 (lzsales1 = gzanywin) lzsales dyrmo*, robust first

#Lagos 2020
data=data.table(read_dta("Lagos2020_data.dta")) #NO X
data=data[(FOMC_Hbased==1|FOMC_Hbased==0)&date>="1994-01-01"&date<="2008-01-01"]
data=data[,c("dr","wr")]
feols(dr~wr,data=data)
#ivregress 2sls ret_m (dr=wr) if (FOMC_Hbased==1|FOMC_Hbased==0)&date>=d(01jan1994)&date<=d(01jan2008), robust
#montests[["Lagos2020"]]=montest(data=data,D="dr",Z="wr",X=NULL,test="simple")
#No X?

##Owens 2020
data=data.table(read_dta("Owens2020_data_1.dta"))
X=c("log_dist_to_park","log_dist_to_highway","log_dist_to_airport","log_dist_to_water","log_dist_to_college")
feols(as.formula(paste0(c("logRj~logAi",X),collapse="+")),data=data)
cols=c(X,"logRj","logAi","logamenities")
data=data[,..cols]
montests[["Owens2020"]]=montest(data=data,X=X,D="logRj",Z="logAi",test="simple")
#ivreg2 logamenities log_dist_to_park log_dist_to_highway log_dist_to_airport log_dist_to_water log_dist_to_college (logRj = logAi), robust

##Pons 2018
data=data.table(read_dta("Pons2018_data.dta"))
data=data[merge_results12==1]
data=data[,c("allocated","treatment","stratum_identifier")]
montests[["Pons2018"]]=montest(data=data,X="stratum_identifier",D="allocated",Z="treatment",test="simple")
#xtivreg2 prop_turnout_pr12t1_an (allocated = treatment) if merge_results12 == 1, i(stratum_identifier) fe robust
