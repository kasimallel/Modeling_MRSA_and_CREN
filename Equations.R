# Load required library
library(ggplot2)
library(deSolve)
library(reshape)
library(dplyr)
library(patchwork)
library(grid) 
library(reshape2)
library(ggpattern)
library(scales) 
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/AMR_BSI_BurdenChile/0_Article_AMR Transmission dynamics Chile & LMICs/0_BSIModelling/0_analysis")

##################################################################################################
#.  MRSA     #
##################################################################################################
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- ---   TIME  --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# Define the time span for simulation
times <- seq(from=0, to=365, by = 1)  # Simulate over a year
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- ---  MODEL PARAMETERS BELOW --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
#Staphylococcus aureus[MRSA/MSSA]####### 
# ----------------------------------------#
zeta3_m_p1 <- 0.231  # mortality rate from IMS  , male
zeta3_f_p1 = zeta3_m_p1*2.07  # mortality rate from IMS  , female
zeta1_m_p1 = zeta3_m_p1*1.01  # mortality rate from IMR, male
zeta1_f_p1 = zeta3_m_p1*1.22  # mortality rate from IMR, female
zeta2_m_p1 = zeta3_m_p1*2.32  # mortality rate from ISR, male
zeta2_f_p1 = zeta3_m_p1*2.25  # mortality rate from ISR, female
zeta4_m_p1 = zeta3_m_p1*1.10  # mortality rate from ISS , male
zeta4_f_p1 = zeta3_m_p1*2.34  # mortality rate from ISS , female
#(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
kappa_p1 = (35.9-6.2)/(49.5-4.85)
parameters1 <- c(
  delta1_p1 = 0.0016,  
  delta2_p1 = 0.0016,  
  Disch_U_f_p1 = 1/6,
  Disch_U_m_p1 = 1/6,
  Disch_CR_f_p1 = 1/6,
  Disch_CR_m_p1 = 1/6,
  Disch_CS_f_p1 = 1/6,
  Disch_CS_m_p1 = 1/6,
  mu0_p1 = 0.52,  # % of women among U
  mu1_p1 = 0.52,  # % of women among CR
  mu2_p1 = 0.52,  # % of women among CS
  mu3_p1 = 0.3165,  # % of women among IMR
  mu4_p1 = 0.3889,  # % of women among ISR
  mu5_p1 = 0.3938,  # % of women among IMS
  mu6_p1 = 0.4487,  # % of women among IMS
  #Exposure to anbiotics; treatment amongst susceptible populations
  psi_m_p1 = 0.1474, # % of  individuals exposed  to vancomycin/penicillin among males
  psi_w_p1 = 0.184,  # % of  individuals exposed to vancomycin/penicillin among females
  #Percentage of people under treatment for MRSA decolonisation
  psi_mtr_p1=0.109, #IEAT is 1.35 times higher among CR; hence if treatment is psi_m_p1; psi_mtr_p1=psi_m_p1/1.35
  psi_wtr_p1=0.136,
  #Fitness cost. c reduces the transmission rate among resistant strains [1/unit time] [%].
  c_p1=0.09,
  #Progression to the development of infection from colonisation among CR and CS states. 
  beta1_m_p1 = (1/21)*0.26,  # inverse of LOS plus progression from colonisation to infection among males CR
  beta2_m_p1 = (1/11)*0.099,  # inverse of LOS plus progression from colonisation to infection among males CS
  beta1_f_p1 = (1/29)*0.26, # inverse of LOS plus progression from colonisation to infection among females CR
  beta2_f_p1 = (1/14)*0.099,  # inverse of LOS plus progression from colonisation to infection among females CS
  #Natural clearance of mild and severe infections among CR and CS states, respectively [1/unit time] [%].
  gamma1_p1 = 0.001,  # Natural clearance among mild infections R
  gamma2_p1 = 0.001,  # Natural clearance among severe infections R
  gamma3_p1 = 0.001,  # Natural clearance among mild infections S
  gamma4_p1 = 0.001, # Natural clearance among severe infections S
  #Mean time of infection considering length of hospital stays [1/length of hospital stay]. 
  omega1_d_m_p1 = (1/20), #IMR patients who died, male
  omega1_r_m_p1 = (1/23), #IMR patients who recovered, male
  omega1_d_f_p1 = (1/13), #IMR patients who died, female
  omega1_r_f_p1 = (1/26), #IMR patients who recovered, female
  omega2_d_m_p1 = (1/11),  #ISR patients who died, male
  omega2_r_m_p1 = (1/18),  #ISR patients who recovered, male
  omega2_d_f_p1 = (1/14), #ISR patients who died, female 
  omega2_r_f_p1 = (1/19), #ISR patients who recovered, female
  omega3_d_m_p1 = (1/12), #IMS patients who died, male
  omega3_r_m_p1 = (1/11), #IMS patients who recovered, male
  omega3_d_f_p1 = (1/20), #IMS patients who died, female
  omega3_r_f_p1 = (1/16), #IMS patients who recovered, female
  omega4_d_m_p1 = (1/11), #ISS patients who died, male
  omega4_r_m_p1 = (1/19), #ISS patients who recovered, male
  omega4_d_f_p1 = (1/14), #ISS patients who died, female
  omega4_r_f_p1 = (1/17), #ISS patients who recovered, female
  #Percentage of inpatients with CR or CS, respectively, progressing to severe infection in intensive care units [1/unit time] [%].
  alpha1_f_p1 = 0.4196, #% patients with CR progressing to severe infection, males
  alpha2_f_p1 = 0.3438, #% patients with CS progressing to severe infection, males
  alpha1_m_p1 = 0.3517, #% patients with CR progressing to severe infection, females
  alpha2_m_p1 = 0.2826, #% patients with CS progressing to severe infection, females
  #Progression from mild to severe infection from IMR and IMS, respectively [1/unit time] [%].
  epsilon1_p1 = 0.01, #progression from IMR to ISR
  epsilon2_p1 = 0.01, #progression from IMS to ISS
  #Mortality rates from infection. ζ1 and ζ2 are mortality rates from mild and severe resistant infections, respectively. ζ3 and ζ4 are from mild and severe susceptible infections, respectively [1/unit time] [%].
  zeta3_m_p1 = 0.231,  # mortality rate from IMS  , male
  zeta3_f_p1 = 0.231*2.07,  # mortality rate from IMS  , female
  zeta1_m_p1 = 0.231*1.01,  # mortality rate from IMR, male
  zeta1_f_p1 = 0.231*1.22,  # mortality rate from IMR, female
  zeta2_m_p1 = 0.231*2.32,  # mortality rate from ISR, male
  zeta2_f_p1 = 0.231*2.25,  # mortality rate from ISR, female
  zeta4_m_p1 = 0.231*1.10,  # mortality rate from ISS , male
  zeta4_f_p1 = 0.231*2.34,  # mortality rate from ISS , female
  #Recovery rates from infection, including IMR, ISR, IMS and ISS due to treatment received [1/unit time] [%].
  nu1_m_p1 = (1-zeta1_m_p1), #recovery rates from IMR, males
  nu1_f_p1 = (1-zeta1_f_p1), #recovery rates from IMR, females
  nu2_m_p1 = (1-zeta2_m_p1), #recovery rates from ISR, males
  nu2_f_p1 = (1-zeta2_f_p1), #recovery rates from ISR, females
  nu3_m_p1 = (1-zeta3_m_p1), #recovery rates from IMS, males
  nu3_f_p1 = (1-zeta3_f_p1), #recovery rates from IMS, females
  nu4_m_p1 =  (1-zeta4_m_p1), #recovery rates from ISS, males
  nu4_f_p1 =  (1-zeta4_f_p1),#recovery rates from ISS, females
  #Constant background rate that captures transmission from non-human sources, horizontal transmission, or de novo emergence [1/unit time] [number].
  b_p1 = 0.01,
  # Percentage of people with resistant infections receiving inappropriate empirical antibiotic treatment [1/unit time] [%].
  phi_m_p1 = 0.459,  # 
  phi_f_p1 = 0.413,  # 
  #Factor of burden associated with inappropriate empirical antibiotic treatment and increased ICU admission among resistant infections [1/unit time] [%].
  pi_p1= 1.35,
  #Transmission parameter {update this correspondingly after calibrating it with real data}
  #tau_p1= 0.03461113, #estimated 
  tau_p1=0.2229124,
  #community-acquired infection upon hospital admission rate
  caIha_p1= 0.001,
  ###
  
  #percentage of people tested
  test_p1=0.20, 
  HR_perc1=0.2,
  or_HR_scenar1_a=1.04,
  or_HR_scenarMen_a=1.37,
  ### ### ### ### ###
  #sensitivity chrom_1
  sens_chrom_a=0.826,
  #sensitivity chrom_1
  sens_chrom2_a=0.622, 
  #sensitivity chrom_1
  sens_pcr_a=0.881,
  #turnaround chrom_1
  turn_chrom_a=3,
  #turnaround chrom_1
  turn_chrom2_a=2,
  #turnaround pcr_1
  turn_pcr_a=1,
  #isolation contact precaution transmission reduction
  reduc_conpre_a=0.365,
  #efficiency decolonisation
  eff_decol_a= 0.53,
  #effect on self-infection decolonisation
  eff_decol_selfi_a=0.32,
  #Turnaround decolonisation program in days
  turnaround_decol_a=5,
  ##
  #costs wards
  c_general_ward= 50,
  c_intermediate_ward=92,
  c_icu_ward=218,
  c_decol_1pd=6.5,
  c_isolation=42.3,
  c_chrom=10.2,
  c_chrom2=13.6, 
  c_pcr=33,
  c_bc=16.9,
  #utilities
  u_healthy=0.92,
  u_icu=0.58,
  u_gw=0.64,
  u_recovICU=0.74
)




# ----------------------------------------#
######
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- ---BASELINE CONDITIONS BELOW--- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
#Staphylococcus aureus[MRSA/MSSA]####### 
# ----------------------------------------#
N <- 1000  # Total population size
# Initial conditions (population sizes in each group)
U_m10 <- 0.7 * N *(1-0.52)
CR_m10 <- 0.1194 * N *(1-0.52)
CS_m10 <-  0.18156* N *(1-0.52)
IMR_m10 <- 0.26 * CR_m10* (1-0.4196)
ISR_m10 <- 0.26 * CR_m10* 0.4196
IMS_m10 <- 0.099 * CS_m10*(1-0.3438)
ISS_m10 <- 0.099 * CS_m10*0.3438
RR_m10 <-0
RS_m10  <-0
DR_m10  <-0
DS_m10 <-0
N_to0 <- 1050
utility_to0<-0
cost_to0<-0
new_admin0<-0
discharge0<- 0

U_f10 <- 0.7 * N *0.52
CR_f10 <- 0.1194 * N *0.52
CS_f10 <-  0.18156* N *0.52
IMR_f10 <- 0.26 * CR_f10* (1-0.3517)
ISR_f10 <- 0.26 * CR_f10*0.3517
IMS_f10 <- 0.099 * CS_f10*(1-0.2826)
ISS_f10 <- 0.099 * CS_f10*0.2826
RR_f10 <-0
RS_f10  <-0
DR_f10  <-0
DS_f10 <-0

N_0m10<-  U_m10 + CR_m10 + CS_m10 + IMR_m10 + ISR_m10 + IMS_m10 + ISS_m10 +  RR_m10 + RS_m10 + DR_m10 +DS_m10
N_0f10<-  U_f10 + CR_f10 + CS_f10 + IMR_f10 + ISR_f10 + IMS_f10 + ISS_f10 +  RR_f10 + RS_f10 + DR_f10 +DS_f10

state1 <- c(U_m1 = U_m10, CR_m1=CR_m10, CS_m1= CS_m10, IMR_m1= IMR_m10, ISR_m1=ISR_m10, IMS_m1= IMS_m10, ISS_m1= ISS_m10, RR_m1= RR_m10, RS_m1=RS_m10, DR_m1= DR_m10, DS_m1=DS_m10,
            U_f1 = U_f10, CR_f1=CR_f10, CS_f1= CS_f10, IMR_f1= IMR_f10, ISR_f1=ISR_f10, IMS_f1= IMS_f10, ISS_f1= ISS_f10, RR_f1= RR_f10, RS_f1=RS_f10, DR_f1= DR_f10, DS_f1=DS_f10, N_to=N_to0, utility=utility_to0, cost=cost_to0, new_admin=new_admin0, discharge= discharge0)

N_orig1<-N_0m10 + N_0f10
N_tdif <- N_orig1

# ----------------------------------#
#####
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- MODEL & EQUATIONS BELOW --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# Define a function for the differential equations
ARB_model_1ch_do_nothing <- function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<- state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2))
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec)+ b_p1*(r_v2)) + ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a+1))
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1)-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+INF_CR_m_p1+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+INF_CS_m_p1+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1)
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1)
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(INF_IMS_m_p1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)+(INF_ISS_m_p1)
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1)-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+INF_CR_f_p1+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) + u_healthy*(RR_f1+RR_m1+RS_f1+RS_m1)      
  dcost <-  c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
#I.1. test+ treatment decolonisation, all new admissions 
ARB_model_1ch_td_newadm <- function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2))
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec)+ b_p1*(r_v2)) + ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a))
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1 + (INF_CR_m_p1*(interv_inf_Rpd1)) +(INF_IMR_m_p1*(interv_inf_Rpd1)) + (INF_ISR_m_p1*(interv_inf_Rpd1)))-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+(INF_CR_m_p1*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_m_p1 +interv_inf_reductPr*INF_ISR_m_p1)+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+(INF_CS_m_p1)+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(INF_IMS_m_p1) 
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)+(INF_ISS_m_p1) 
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1+ (INF_CR_f_p1*(interv_inf_Rpd1)) +(INF_IMR_f_p1*(interv_inf_Rpd1)) + (INF_ISR_f_p1*(interv_inf_Rpd1)) )-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+(INF_CR_f_p1*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_f_p1 +interv_inf_reductPr*INF_ISR_f_p1)+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) + u_healthy*(RR_f1+RR_m1+RS_f1+RS_m1)      
  dcost <- (influx_nonARB)*(c_chrom) + (influx_ARB)*(c_chrom+c_decol_1pd) + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<- influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
ARB_model_1ch2_td_newadm <-function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2))
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec)+ b_p1*(r_v2)) + ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  interv_inf_Rpd1<- sens_chrom2_a*eff_decol_a*(1/(turn_chrom2_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom2_a))
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1 + (INF_CR_m_p1*(interv_inf_Rpd1)) +(INF_IMR_m_p1*(interv_inf_Rpd1)) + (INF_ISR_m_p1*(interv_inf_Rpd1)))-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+(INF_CR_m_p1*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_m_p1 +interv_inf_reductPr*INF_ISR_m_p1 )+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+(INF_CS_m_p1)+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(INF_IMS_m_p1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)+(INF_ISS_m_p1) 
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1+ (INF_CR_f_p1*(interv_inf_Rpd1)) +(INF_IMR_f_p1*(interv_inf_Rpd1)) + (INF_ISR_f_p1*(interv_inf_Rpd1)))-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+(INF_CR_f_p1*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_f_p1 +interv_inf_reductPr*INF_ISR_f_p1 +interv_inf_reductPr)+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (influx_nonARB)*(c_chrom2) + (influx_ARB)*(c_chrom2+c_decol_1pd) + c_general_ward*(U_m1+ CR_m1+ CS_m1+ U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<- influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
ARB_model_1pcr_td_newadm <-function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  discharge<-state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"] 
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2))
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec)+ b_p1*(r_v2)) + ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  interv_inf_Rpd1<- sens_pcr_a*eff_decol_a*(1/(turn_pcr_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_pcr_a))
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1 + (INF_CR_m_p1*(interv_inf_Rpd1)) +(INF_IMR_m_p1*(interv_inf_Rpd1)) + (INF_ISR_m_p1*(interv_inf_Rpd1)))-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+(INF_CR_m_p1*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_m_p1 +interv_inf_reductPr*INF_ISR_m_p1 )+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+(INF_CS_m_p1)+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(INF_IMS_m_p1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)+(INF_ISS_m_p1) 
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1+ (INF_CR_f_p1*(interv_inf_Rpd1)) +(INF_IMR_f_p1*(interv_inf_Rpd1)) + (INF_ISR_f_p1*(interv_inf_Rpd1)))-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+(INF_CR_f_p1*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_f_p1 +interv_inf_reductPr*INF_ISR_f_p1 +interv_inf_reductPr)+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (influx_nonARB)*(c_pcr) + (influx_ARB)*(c_pcr+c_decol_1pd) + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<-influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin,ddischarge))
  
  return(results1)
  
}
#I.2. test+treatment isolation, all new admissions 
ARB_model_1ch_tiso_newadm <-function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<- state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  #interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  #interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a))
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*(1-(sens_chrom_a*(1/(turn_chrom_a))*reduc_conpre_a))*(((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1))))/Nt1_spec) + b_p1*(r_v2)) 
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  FOC_cr_1 + FOC_cs_1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1)-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+INF_CR_m_p1+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+INF_CS_m_p1+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1)
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1)
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(INF_IMS_m_p1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)+(INF_ISS_m_p1)
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1)-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+INF_CR_f_p1+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (influx_nonARB)*(c_chrom) + (influx_ARB)*(c_chrom+ c_isolation)  + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<-influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
ARB_model_1ch2_tiso_newadm <-function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  ddischarge<- state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  #interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  #interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a+1))
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*(1-(sens_chrom2_a*(1/(turn_chrom2_a))*reduc_conpre_a))*(((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1))))/Nt1_spec) + b_p1*(r_v2)) 
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  FOC_cr_1 + FOC_cs_1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1)-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+INF_CR_m_p1+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+INF_CS_m_p1+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1)
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1)
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(INF_IMS_m_p1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)+(INF_ISS_m_p1)
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1)-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+INF_CR_f_p1+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (influx_nonARB)*(c_chrom2) + (influx_ARB)*(c_chrom2+ c_isolation) + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<-influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
} 
ARB_model_1pcr_tiso_newadm <- function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<- state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  #interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  #interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a+1))
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  
  # DEFINITION OF THE FORCE OF INFECTION
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*(1-(sens_pcr_a*(1/(turn_pcr_a))*reduc_conpre_a))*(((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1))))/Nt1_spec) + b_p1*(r_v2)) 
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  FOC_cr_1 + FOC_cs_1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1)-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+INF_CR_m_p1+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+INF_CS_m_p1+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1)
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1)
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(INF_IMS_m_p1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)+(INF_ISS_m_p1)
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1)-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+INF_CR_f_p1+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (influx_nonARB)*(c_pcr) + (influx_ARB)*(c_pcr+ c_isolation)  + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<-influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
#II.1. treatment decolonisation, high-risk patients: males
ARB_model_1ch_td_newadmHR_m <- function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  new_admin<-state["new_admin"]
  discharge<- state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  HR_perc1<- parms["HR_perc1"]
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2))
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec)+ b_p1*(r_v2)) + ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a)) #clearance
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a)) #infection self reduction
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  Influx_men <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1 +INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_r <-  INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_nor <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1
  Influx_women <-INF_U_f_p1 + INF_CS_f_p1+ INF_IMS_f_p1 + INF_ISS_f_p1+INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_r <-  INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_nor <-  INF_U_f_p1  + INF_CS_f_p1 + INF_IMS_f_p1 + INF_ISS_f_p1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1 + (INF_CR_m_p1*(interv_inf_Rpd1)) +(INF_IMR_m_p1*(interv_inf_Rpd1)) + (INF_ISR_m_p1*(interv_inf_Rpd1)))-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+(INF_CR_m_p1*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_m_p1 +interv_inf_reductPr*INF_ISR_m_p1)+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+(INF_CS_m_p1)+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1)-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+(INF_CR_f_p1)+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (Influx_men_nor)*(c_chrom) + (Influx_men_r)*(c_chrom+c_decol_1pd) + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<-influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
ARB_model_1ch2_td_newadmHR_m <- function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  HR_perc1<- parms["HR_perc1"]
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2))
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec)+ b_p1*(r_v2)) + ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  interv_inf_Rpd1<- sens_chrom2_a*eff_decol_a*(1/(turn_chrom2_a+turnaround_decol_a)) #clearance
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom2_a)) #infection self reduction
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  Influx_men <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1 +INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_r <-  INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_nor <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1
  Influx_women <-INF_U_f_p1 + INF_CS_f_p1+ INF_IMS_f_p1 + INF_ISS_f_p1+INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_r <-  INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_nor <-  INF_U_f_p1  + INF_CS_f_p1 + INF_IMS_f_p1 + INF_ISS_f_p1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1 + (INF_CR_m_p1*(interv_inf_Rpd1)) +(INF_IMR_m_p1*(interv_inf_Rpd1)) + (INF_ISR_m_p1*(interv_inf_Rpd1)))-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+(INF_CR_m_p1*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_m_p1 +interv_inf_reductPr*INF_ISR_m_p1)+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+(INF_CS_m_p1)+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1)-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+(INF_CR_f_p1)+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (Influx_men_nor)*(c_chrom2) + (Influx_men_r)*(c_chrom2+c_decol_1pd) + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<-influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
ARB_model_1_pcr_td_newadmHR_m <- function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  HR_perc1<- parms["HR_perc1"]
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2))
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec)+ b_p1*(r_v2)) + ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  interv_inf_Rpd1<- sens_pcr_a*eff_decol_a*(1/(turn_pcr_a+turnaround_decol_a)) #clearance
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_pcr_a)) #infection self reduction
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  Influx_men <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1 +INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_r <-  INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_nor <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1
  Influx_women <-INF_U_f_p1 + INF_CS_f_p1+ INF_IMS_f_p1 + INF_ISS_f_p1+INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_r <-  INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_nor <-  INF_U_f_p1  + INF_CS_f_p1 + INF_IMS_f_p1 + INF_ISS_f_p1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1 + (INF_CR_m_p1*(interv_inf_Rpd1)) +(INF_IMR_m_p1*(interv_inf_Rpd1)) + (INF_ISR_m_p1*(interv_inf_Rpd1)))-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+(INF_CR_m_p1*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_m_p1 +interv_inf_reductPr*INF_ISR_m_p1)+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+(INF_CS_m_p1)+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1)-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+(INF_CR_f_p1)+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (Influx_men_nor)*(c_pcr) + (Influx_men_r)*(c_pcr+c_decol_1pd) + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<-influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
#II.2. treatment decolonisation, high-risk patients: females
ARB_model_1ch_td_newadmHR_f <- function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  HR_perc1<- parms["HR_perc1"]
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2))
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec)+ b_p1*(r_v2)) + ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a)) #clearance
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a)) #infection self reduction
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  Influx_men <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1 +INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_r <-  INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_nor <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1
  Influx_women <-INF_U_f_p1 + INF_CS_f_p1+ INF_IMS_f_p1 + INF_ISS_f_p1+INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_r <-  INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_nor <-  INF_U_f_p1  + INF_CS_f_p1 + INF_IMS_f_p1 + INF_ISS_f_p1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1)-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+(INF_CR_m_p1)+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+(INF_CS_m_p1)+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1)
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1)
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(INF_IMS_m_p1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)+(INF_ISS_m_p1)
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1 + (INF_CR_f_p1*(interv_inf_Rpd1)) +(INF_IMR_f_p1*(interv_inf_Rpd1)) + (INF_ISR_f_p1*(interv_inf_Rpd1)))-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)+ (INF_CR_f_p1*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_f_p1 +interv_inf_reductPr*INF_ISR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (Influx_women_nor)*(c_chrom) + (Influx_women_r)*(c_chrom+c_decol_1pd) + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<-influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
ARB_model_1ch2_td_newadmHR_f <- function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  new_admin<-state["new_admin"]
  discharge<- state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  HR_perc1<- parms["HR_perc1"]
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2))
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec)+ b_p1*(r_v2)) + ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  interv_inf_Rpd1<- sens_chrom2_a*eff_decol_a*(1/(turn_chrom2_a+turnaround_decol_a)) #clearance
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom2_a)) #infection self reduction
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  Influx_men <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1 +INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_r <-  INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_nor <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1
  Influx_women <-INF_U_f_p1 + INF_CS_f_p1+ INF_IMS_f_p1 + INF_ISS_f_p1+INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_r <-  INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_nor <-  INF_U_f_p1  + INF_CS_f_p1 + INF_IMS_f_p1 + INF_ISS_f_p1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1)-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+(INF_CR_m_p1)+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+(INF_CS_m_p1)+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1)
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1)
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(INF_IMS_m_p1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)+(INF_ISS_m_p1)
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1 + (INF_CR_f_p1*(interv_inf_Rpd1)) +(INF_IMR_f_p1*(interv_inf_Rpd1)) + (INF_ISR_f_p1*(interv_inf_Rpd1)))-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)+ (INF_CR_f_p1*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_f_p1 +interv_inf_reductPr*INF_ISR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (Influx_women_nor)*(c_chrom2) + (Influx_women_r)*(c_chrom2+c_decol_1pd) + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<-influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
ARB_model_1_pcr_td_newadmHR_f <- function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  HR_perc1<- parms["HR_perc1"]
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2))
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  (((tau_p1*(1-c_p1)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec)+ b_p1*(r_v2)) + ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  interv_inf_Rpd1<- sens_pcr_a*eff_decol_a*(1/(turn_pcr_a+turnaround_decol_a)) #clearance
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_pcr_a)) #infection self reduction
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  Influx_men <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1 +INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_r <-  INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_nor <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1
  Influx_women <-INF_U_f_p1 + INF_CS_f_p1+ INF_IMS_f_p1 + INF_ISS_f_p1+INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_r <-  INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_nor <-  INF_U_f_p1  + INF_CS_f_p1 + INF_IMS_f_p1 + INF_ISS_f_p1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1)-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+(INF_CR_m_p1)+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+(INF_CS_m_p1)+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1)
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1)
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(INF_IMS_m_p1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)+(INF_ISS_m_p1)
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1 + (INF_CR_f_p1*(interv_inf_Rpd1)) +(INF_IMR_f_p1*(interv_inf_Rpd1)) + (INF_ISR_f_p1*(interv_inf_Rpd1)))-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)+ (INF_CR_f_p1*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_f_p1 +interv_inf_reductPr*INF_ISR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (Influx_women_nor)*(c_pcr) + (Influx_women_r)*(c_pcr+c_decol_1pd) + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<-influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
#III.1 Pre-emptive isolation of all new admissions [NO TEST]
ARB_model_1preE_newadm_all <-function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  #interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  #interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a+1))
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  Influx_men <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1 +INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_r <-  INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_nor <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1
  Influx_women <-INF_U_f_p1 + INF_CS_f_p1+ INF_IMS_f_p1 + INF_ISS_f_p1+INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_r <-  INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_nor <-  INF_U_f_p1  + INF_CS_f_p1 + INF_IMS_f_p1 + INF_ISS_f_p1
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  
  #RESCALING TRANSMISSION PARAMETER
  tau_p1_rs <- ((57.312+62.088)*tau_p1) /((1-0.6069)*(62.088)+0.6069*(57.312)) #(CR_m+CR_f1)*tau_p1 /((1-0.6069)*(CR_f1)+0.6069*(CR_m))
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions  0.3*Beta*women*uncolonised/(Nt) + 0.7*(1-Clevel)*Beta*(men)*uncolonised/Nt  Clevel= Coverage*efficacy*OR
  #FOC_cr_1 <- ((((1-0.6069)*tau_p1_rs*(1-c_p1)*((CR_f1+IMR_f1+ISR_f1))*(U_f1+U_m1)))/Nt1_spec)+(((0.6069*tau_p1_rs*(1-c_p1)*((CR_m1+IMR_m1+ISR_m1)*(1-reduc_conpre_a*or_HR_scenarMen_a))*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2)
  #FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  #FOC_u_1 <-  ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*(1-reduc_conpre_a)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2))
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <- FOC_cr_1 + FOC_cs_1
  
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1)-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+INF_CR_m_p1+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+INF_CS_m_p1+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1)
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1)
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(INF_IMS_m_p1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)+(INF_ISS_m_p1)
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1)-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+INF_CR_f_p1+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (influx_nonARB + influx_ARB)*(c_isolation) + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<-influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
#III.2 Pre-emptive isolation of males
ARB_model_1preE_newadm_m <- function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  #interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  #interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a+1))
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  Influx_men <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1 +INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_r <-  INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_nor <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1
  Influx_women <-INF_U_f_p1 + INF_CS_f_p1+ INF_IMS_f_p1 + INF_ISS_f_p1+INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_r <-  INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_nor <-  INF_U_f_p1  + INF_CS_f_p1 + INF_IMS_f_p1 + INF_ISS_f_p1
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  
  #RESCALING TRANSMISSION PARAMETER
  tau_p1_rs <- (57.312+62.088)*tau_p1 /((1-0.6069)*(62.088)+(0.6069*(57.312))) #(CR_m+CR_f1)*tau_p1 /((1-0.6069)*(CR_f1)+0.6069*(CR_m))
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions  0.3*Beta*women*uncolonised/(Nt) + 0.7*(1-Clevel)*Beta*(men)*uncolonised/Nt  Clevel= Coverage*efficacy*OR
  FOC_cr_1 <- ((((1-0.6069)*tau_p1_rs*(1-c_p1)*((CR_f1+IMR_f1+ISR_f1))*(U_f1+U_m1)))/Nt1_spec)+(((0.6069*tau_p1_rs*(1-c_p1)*((CR_m1+IMR_m1+ISR_m1)*(1-reduc_conpre_a*or_HR_scenarMen_a))*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2)
  FOC_cs_1 <- ((tau_p1*((CS_m1+IMS_m1+ISS_m1)*(U_f1+U_m1)))/Nt1_spec)  + ((tau_p1*((CS_f1+IMS_f1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  FOC_cs_1+ FOC_cr_1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1)-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+INF_CR_m_p1+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+INF_CS_m_p1+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1)
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1)
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(INF_IMS_m_p1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)+(INF_ISS_m_p1)
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1)-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+INF_CR_f_p1+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (Influx_men)*(c_isolation) + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<-influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
#III.3 Pre-emptive isolation of females [NO TEST/ no OR]
ARB_model_1preE_newadm_f <- function(times, state, parms){
  
  # Men 
  U_m1 <- state["U_m1"]
  CR_m1 <- state["CR_m1"]
  CS_m1 <- state["CS_m1"]
  IMR_m1 <- state["IMR_m1"]
  ISR_m1 <- state["ISR_m1"]
  IMS_m1 <- state["IMS_m1"]
  ISS_m1 <- state["ISS_m1"]
  RR_m1 <- state["RR_m1"]
  RS_m1 <- state["RS_m1"]
  DR_m1 <- state["DR_m1"]
  DS_m1 <- state["DS_m1"]
  
  N1_1 <- U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1 + RR_m1 + RS_m1 + DR_m1 + DS_m1
  
  # Women   
  U_f1 <- state["U_f1"]
  CR_f1 <- state["CR_f1"]
  CS_f1 <- state["CS_f1"]
  IMR_f1 <- state["IMR_f1"]
  ISR_f1 <- state["ISR_f1"]
  IMS_f1 <- state["IMS_f1"]
  ISS_f1 <- state["ISS_f1"]
  RR_f1 <- state["RR_f1"]
  RS_f1 <- state["RS_f1"]
  DR_f1 <- state["DR_f1"]
  DS_f1 <- state["DS_f1"]
  N_to <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_1 <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 + RR_f1 + RS_f1 + DR_f1 + DS_f1
  
  #N total (women+men)
  Nt_1 <- max(N1_1 + N2_1, 1)
  #population at time t
  Nt1_spec <- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  # # # # # # #
  
  #Extract parameters
  delta1_p1<- parms["delta1_p1"]
  delta2_p1<- parms["delta2_p1"]
  Disch_U_f_p1<-parms["Disch_U_f_p1"]
  Disch_U_m_p1<- parms["Disch_U_m_p1"] 
  Disch_CR_f_p1<- parms["Disch_CR_f_p1"] 
  Disch_CR_m_p1<-parms["Disch_CR_m_p1"] 
  Disch_CS_f_p1<-parms["Disch_CS_f_p1"] 
  Disch_CS_m_p1 <-parms["Disch_CS_m_p1"]
  mu0_p1<- parms["mu0_p1"]
  mu1_p1<- parms["mu1_p1"]
  mu2_p1<- parms["mu2_p1"]
  mu3_p1<- parms["mu3_p1"]
  mu4_p1<- parms["mu4_p1"]
  mu5_p1<- parms["mu5_p1"]
  mu6_p1<- parms["mu6_p1"]
  psi_m_p1<- parms["psi_m_p1"]
  psi_w_p1<- parms["psi_w_p1"]
  c_p1<- parms["c_p1"]
  beta1_m_p1<- parms["beta1_m_p1"]
  beta2_m_p1 <- parms["beta2_m_p1"]
  beta1_f_p1<- parms["beta1_f_p1"]
  beta2_f_p1<- parms["beta2_f_p1"]
  gamma1_p1<- parms["gamma1_p1"]
  gamma2_p1<- parms["gamma2_p1"]
  gamma3_p1<- parms["gamma3_p1"]
  gamma4_p1<- parms["gamma4_p1"]
  omega1_d_m_p1<- parms["omega1_d_m_p1"]
  omega1_r_m_p1<- parms["omega1_r_m_p1"]
  omega1_d_f_p1<- parms["omega1_d_f_p1"]
  omega1_r_f_p1<- parms["omega1_r_f_p1"]
  omega2_d_m_p1<- parms["omega2_d_m_p1"]
  omega2_r_m_p1<- parms["omega2_r_m_p1"]
  omega2_d_f_p1<- parms["omega2_d_f_p1"]
  omega2_r_f_p1<- parms["omega2_r_f_p1"]
  omega3_d_m_p1<- parms["omega3_d_m_p1"]
  omega3_r_m_p1<- parms["omega3_r_m_p1"]
  omega3_d_f_p1<- parms["omega3_d_f_p1"]
  omega3_r_f_p1<- parms["omega3_r_f_p1"]
  omega4_d_m_p1<- parms["omega4_d_m_p1"]
  omega4_r_m_p1<- parms["omega4_r_m_p1"]
  omega4_d_f_p1<- parms["omega4_d_f_p1"]
  omega4_r_f_p1<- parms["omega4_r_f_p1"]
  alpha1_m_p1<- parms["alpha1_m_p1"]
  alpha2_m_p1<- parms["alpha2_m_p1"]
  alpha1_f_p1<- parms["alpha1_f_p1"]
  alpha2_f_p1<- parms["alpha2_f_p1"]
  epsilon1_p1<- parms["epsilon1_p1"]
  epsilon2_p1<- parms["epsilon2_p1"]
  zeta3_m_p1<- parms["zeta3_m_p1"]
  zeta3_f_p1<- parms["zeta3_f_p1"]
  zeta1_m_p1<- parms["zeta1_m_p1"]
  zeta1_f_p1<- parms["zeta1_f_p1"]
  zeta2_m_p1<- parms["zeta2_m_p1"]
  zeta2_f_p1<- parms["zeta2_f_p1"]
  zeta4_m_p1<- parms["zeta4_m_p1"]
  zeta4_f_p1<- parms["zeta4_f_p1"]
  nu1_m_p1<- parms["nu1_m_p1"]
  nu1_f_p1<- parms["nu1_f_p1"]
  nu2_m_p1 <- parms["nu2_m_p1"]
  nu2_f_p1<- parms["nu2_f_p1"]
  nu3_m_p1<- parms["nu3_m_p1"]
  nu3_f_p1<- parms["nu3_f_p1"]
  nu4_m_p1<- parms["nu4_m_p1"]
  nu4_f_p1<- parms["nu4_f_p1"]
  b_p1<- parms["b_p1"]
  phi_m_p1<- parms["phi_m_p1"]
  phi_f_p1<- parms["phi_f_p1"]
  pi_p1<- parms["pi_p1"]
  tau_p1<-parms["tau_p1"]
  caIha_p1<-parms["caIha_p1"]
  psi_mtr_p1 <-parms["psi_mtr_p1"]
  psi_wtr_p1 <-parms["psi_wtr_p1"]
  #percentage of people tested
  test_p1<-parms["test_p1"] 
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #Prevalence of MRSA
  P1_t1 <- (CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1)/(CR_f1 + IMR_f1 + ISR_f1 + CR_m1 + IMR_m1 + ISR_m1 + CS_f1 + IMS_f1 + ISS_f1 + CS_m1 + IMS_m1 + ISS_m1)
  
  #Random value for competing transmissions
  ra_v <- runif(1, min = 0.00, max = 0.01)
  ra_v=0
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v2 <- rbeta(1, alpha12, beta12)
  h_ieat1_p1 <- (alpha1_m_p1)/((pi_p1*phi_m_p1)+(1-phi_m_p1))
  h_ieat2_p1 <- (alpha1_f_p1)/((pi_p1*phi_f_p1)+(1-phi_f_p1))
  N_to<- U_f1 + CR_f1 + CS_f1 + IMR_f1 + ISR_f1 + IMS_f1 + ISS_f1 +U_m1 + CR_m1 + CS_m1 + IMR_m1 + ISR_m1 + IMS_m1 + ISS_m1
  
  #INTERVENTION ADJUSTMENTS:
  #calculation of clearance per day among influx to the hospital being ARB and receiving decol treatment
  #interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  #interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a+1))
  
  #Influx of populations 
  INF_U_f_p1 <- (1050- Nt1_spec)*0.7*mu0_p1
  INF_U_m_p1 <- (1050- Nt1_spec)*0.7*(1-mu0_p1)
  INF_CR_f_p1 <- (1050- Nt1_spec)*0.1194*mu1_p1
  INF_CR_m_p1<- (1050- Nt1_spec)*0.1194*(1-mu1_p1)
  INF_CS_f_p1<- (1050- Nt1_spec)*0.18156*mu2_p1
  INF_CS_m_p1 <- (1050- Nt1_spec)*0.18156*(1-mu2_p1)
  INF_IMR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu3_p1*(1/8)
  INF_IMR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu3_p1)*(1/8)
  INF_ISR_f_p1<- (1050- Nt1_spec)*(caIha_p1)*mu4_p1*(1/8)
  INF_ISR_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu4_p1)*(1/8)
  INF_IMS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu5_p1)*(1/8)
  INF_IMS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu5_p1)*(1/8)
  INF_ISS_f_p1<- (1050- Nt1_spec)*(caIha_p1)*(mu6_p1)*(1/8)
  INF_ISS_m_p1<- (1050- Nt1_spec)*(caIha_p1)*(1-mu6_p1)*(1/8)
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  Influx_men <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1 +INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_r <-  INF_CR_m_p1+INF_IMR_m_p1 +INF_ISR_m_p1
  Influx_men_nor <-  INF_U_m_p1  + INF_CS_m_p1 + INF_IMS_m_p1 + INF_ISS_m_p1
  Influx_women <-INF_U_f_p1 + INF_CS_f_p1+ INF_IMS_f_p1 + INF_ISS_f_p1+INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_r <-  INF_CR_f_p1 +INF_IMR_f_p1 +INF_ISR_f_p1
  Influx_women_nor <-  INF_U_f_p1  + INF_CS_f_p1 + INF_IMS_f_p1 + INF_ISS_f_p1
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p1 + INF_U_m_p1 + INF_CS_f_p1 + INF_CS_m_p1 + INF_IMS_f_p1 + INF_IMS_m_p1 + INF_ISS_f_p1 + INF_ISS_m_p1
  influx_ARB<- INF_CR_f_p1 + INF_CR_m_p1+  INF_IMR_f_p1 + INF_IMR_m_p1 + INF_ISR_f_p1 + INF_ISR_m_p1
  
  #RESCALING TRANSMISSION PARAMETER
  tau_p1_rs <- (57.312+62.088)*tau_p1 /((1-0.6069)*(62.088)+0.6069*(57.312)) #(CR_m+CR_f1)*tau_p1 /((1-0.6069)*(CR_f1)+0.6069*(CR_m))
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions  0.3*Beta*women*uncolonised/(Nt) + 0.7*(1-Clevel)*Beta*(men)*uncolonised/Nt  Clevel= Coverage*efficacy*OR
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*(1-reduc_conpre_a)*((CR_f1+IMR_f1+ISR_f1))*(U_f1+U_m1)))/Nt1_spec)+(((tau_p1*(1-c_p1)*((CR_m1+IMR_m1+ISR_m1))*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2)
  FOC_cs_1 <- ((tau_p1*((CS_f1+IMS_f1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec) + ((tau_p1*((CS_f1+IMS_f1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  FOC_u_1 <-  FOC_cs_1 + FOC_cr_1
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m1 <-   (delta1_p1*CR_m1)+(delta2_p1*CS_m1)+(INF_U_m_p1)-(FOC_u_1*(1-mu0_p1))-(U_m1*Disch_U_m_p1)+(psi_m_p1*CS_m1)+(psi_mtr_p1*CR_m1)
  dCR_m1 <- -(delta1_p1*CR_m1)-(beta1_m_p1*CR_m1)-(psi_mtr_p1*CR_m1)+(gamma1_p1*IMR_m1)+(gamma2_p1*ISR_m1)+INF_CR_m_p1+((1-mu1_p1)*(FOC_cr_1))-(CR_m1*Disch_CR_m_p1)
  dCS_m1 <- -(delta2_p1*CS_m1)-(beta2_m_p1*CS_m1)-(psi_m_p1*CS_m1)  +(gamma3_p1*IMS_m1)+(gamma4_p1*ISS_m1)+INF_CS_m_p1+((1-mu2_p1)*(FOC_cs_1))-(CS_m1*Disch_CS_m_p1)
  dIMR_m1 <- ((beta1_m_p1*CR_m1)*(1-alpha1_m_p1))-(gamma1_p1*IMR_m1)-(omega1_r_m_p1*nu1_m_p1*IMR_m1)-(epsilon1_p1*IMR_m1)-(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(INF_IMR_m_p1)
  dISR_m1 <- (beta1_m_p1*CR_m1*alpha1_m_p1)      -(gamma2_p1*ISR_m1)-(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(epsilon1_p1*IMR_m1)-(omega2_d_m_p1*zeta2_m_p1*ISR_m1)+(INF_ISR_m_p1)
  dIMS_m1 <- (beta2_m_p1*CS_m1*(1-alpha2_m_p1))  -(gamma3_p1*IMS_m1)-(omega3_r_m_p1*nu3_m_p1*IMS_m1)-(epsilon2_p1*ISS_m1)-(omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(INF_IMS_m_p1)
  dISS_m1 <- (beta2_m_p1*CS_m1*(alpha2_m_p1))    -(gamma4_p1*ISS_m1)-(omega4_r_m_p1*nu4_m_p1*ISS_m1)+(epsilon2_p1*ISS_m1)-(omega4_d_m_p1*zeta4_m_p1*ISS_m1)+(INF_ISS_m_p1)
  dRR_m1 <-  (omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)
  dRS_m1 <-  (omega3_r_m_p1*nu3_m_p1*IMS_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)
  dDR_m1 <-  (omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  dDS_m1 <-  (omega3_d_m_p1*zeta3_m_p1*IMS_m1)+(omega4_d_m_p1*zeta4_m_p1*ISS_m1)
  
  dU_f1 <-   (delta1_p1*CR_f1)+(delta2_p1*CS_f1)+(INF_U_f_p1)-(FOC_u_1*mu0_p1)-(U_f1*Disch_U_f_p1)+(psi_w_p1*CS_f1)+(psi_wtr_p1*CR_f1)
  dCR_f1 <- -(delta1_p1*CR_f1)-(beta1_f_p1*CR_f1)-(psi_wtr_p1*CS_f1)+(gamma1_p1*IMR_f1)+(gamma2_p1*ISR_f1)+INF_CR_f_p1+((mu1_p1)*(FOC_cr_1))-(CR_f1*Disch_CR_f_p1)
  dCS_f1<-  -(delta2_p1*CS_f1)-(beta2_f_p1*CS_f1)-(psi_w_p1*CS_f1)  +(gamma3_p1*IMS_f1)+(gamma4_p1*ISS_f1)+INF_CS_f_p1+((mu2_p1)*(FOC_cs_1))-(CS_f1*Disch_CS_f_p1)
  dIMR_f1 <- ((beta1_f_p1*CR_f1)*(1-alpha1_f_p1))-(gamma1_p1*IMR_f1)-(omega1_r_f_p1*nu1_f_p1*IMR_f1)-(epsilon1_p1*IMR_f1)-(omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(INF_IMR_f_p1)
  dISR_f1 <- (beta1_f_p1*CR_f1*alpha1_f_p1)      -(gamma2_p1*ISR_f1)-(omega2_r_f_p1*nu2_f_p1*ISR_f1)+(epsilon1_p1*IMR_f1)-(omega2_d_f_p1*zeta2_f_p1*ISR_f1)+(INF_ISR_f_p1)
  dIMS_f1 <- (beta2_f_p1*CS_f1*(1-alpha2_f_p1))  -(gamma3_p1*IMS_f1)-(omega3_r_f_p1*nu3_f_p1*IMS_f1)-(epsilon2_p1*ISS_f1)-(omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(INF_IMS_f_p1)
  dISS_f1 <- (beta2_f_p1*CS_f1*(alpha2_f_p1))    -(gamma4_p1*ISS_f1)-(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(epsilon2_p1*ISS_f1)-(omega4_d_f_p1*zeta4_f_p1*ISS_f1)+(INF_ISS_f_p1)
  dRR_f1 <- (omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega2_r_f_p1*nu2_f_p1*ISR_f1)
  dRS_f1 <- (omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)
  dDR_f1 <- (omega1_d_f_p1*zeta1_f_p1*IMR_f1)+(omega2_d_f_p1*zeta2_f_p1*ISR_f1)
  dDS_f1 <- (omega3_d_f_p1*zeta3_f_p1*IMS_f1)+(omega4_d_f_p1*zeta4_f_p1*ISS_f1)
  dN_to<- dU_m1+ dCR_m1+ dCS_m1+ dIMR_m1+ dISR_m1+ dIMS_m1+ dISS_m1 +dU_f1+ dCR_f1+ dCS_f1+ dIMR_f1+ dISR_f1+ dIMS_f1+ dISS_f1
  dutility <- u_healthy*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +u_icu*(ISR_m1+ISS_m1+ISR_f1+ISS_f1)+ u_gw*(IMS_f1+ IMR_f1+IMS_m1+ IMR_m1) +u_recovICU*((omega2_r_f_p1*nu2_f_p1*ISR_f1)+(omega4_r_f_p1*nu4_f_p1*ISS_f1)+(omega2_r_m_p1*nu2_m_p1*ISR_m1)+(omega4_r_m_p1*nu4_m_p1*ISS_m1)) + u_healthy*((omega1_r_f_p1*nu1_f_p1*IMR_f1)+(omega3_r_f_p1*nu3_f_p1*IMS_f1)+(omega1_r_m_p1*nu1_m_p1*IMR_m1)+(omega3_r_m_p1*nu3_m_p1*IMS_m1))      
  dcost <- (Influx_women)*(c_isolation) + c_general_ward*(U_m1+ CR_m1+ CS_m1+U_f1+ CR_f1+ CS_f1) +c_intermediate_ward*(IMR_m1+ IMS_m1+IMR_f1+ IMS_f1)+ c_icu_ward*(ISR_m1+ ISS_m1+ISR_f1+ ISS_f1) 
  dnew_admin<-influx_nonARB + influx_ARB
  ddischarge <- U_m1*Disch_U_m_p1+CR_m1*Disch_CR_m_p1+CS_m1*Disch_CS_m_p1+U_f1*Disch_U_f_p1+CR_f1*Disch_CR_f_p1+CS_f1*Disch_CS_f_p1
  #discharge<- state["discharge"] #list results ddischarge
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to, dutility, dcost, dnew_admin, ddischarge))
  
  return(results1)
  
}
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- SOLVING EQUATIONS--- ------ --- --- ------ --- --- --- --- --- #
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
#TEST_odes#######
#ode(y = state1, times = times, func = ARB_model_1ch_do_nothing, parms = parameters1, method = "rk4")
#ode(y = state1, times = times, func = ARB_model_1ch_td_newadm, parms = parameters1, method = "rk4")
#ode(y = state1, times = times, func = ARB_model_1ch2_td_newadm, parms = parameters1, method = "rk4")
#ode(y = state1, times = times, func = ARB_model_1pcr_td_newadm, parms = parameters1, method = "rk4")

#ode(y = state1, times = times, func = ARB_model_1ch_tiso_newadm, parms = parameters1, method = "rk4")
#ode(y = state1, times = times, func = ARB_model_1ch2_tiso_newadm, parms = parameters1, method = "rk4")
#ode(y = state1, times = times, func = ARB_model_1pcr_tiso_newadm, parms = parameters1, method = "rk4")

#ode(y = state1, times = times, func = ARB_model_1ch_td_newadmHR_m, parms = parameters1, method = "rk4")
#ode(y = state1, times = times, func = ARB_model_1ch2_td_newadmHR_m, parms = parameters1, method = "rk4")
#ode(y = state1, times = times, func = ARB_model_1_pcr_td_newadmHR_m, parms = parameters1, method = "rk4")

#ode(y = state1, times = times, func = ARB_model_1ch_td_newadmHR_f, parms = parameters1, method = "rk4")
#ode(y = state1, times = times, func = ARB_model_1ch2_td_newadmHR_f, parms = parameters1, method = "rk4")
#ode(y = state1, times = times, func = ARB_model_1_pcr_td_newadmHR_f, parms = parameters1, method = "rk4")

#O_solution<-ode(y = state1, times = times, func = ARB_model_1preE_newadm_all , parms = parameters1, method = "rk4")
#O_solution <-ode(y = state1, times = times, func = ARB_model_1preE_newadm_m, parms = parameters1, method = "rk4")
#O_solution <-ode(y = state1, times = times, func = ARB_model_1preE_newadm_f, parms = parameters1, method = "rk4")

#O_solution<- as.data.frame(O_solution)
#O_solution[366,"CR_f1"]
#####

# Define the models you want to run according to the functions (per strategy)
model_list <- c("ARB_model_1ch_do_nothing", "ARB_model_1ch_td_newadm","ARB_model_1ch2_td_newadm","ARB_model_1pcr_td_newadm","ARB_model_1ch_tiso_newadm","ARB_model_1ch2_tiso_newadm","ARB_model_1pcr_tiso_newadm",
                "ARB_model_1ch_td_newadmHR_m","ARB_model_1ch2_td_newadmHR_m","ARB_model_1_pcr_td_newadmHR_m","ARB_model_1ch_td_newadmHR_f","ARB_model_1ch2_td_newadmHR_f","ARB_model_1_pcr_td_newadmHR_f","ARB_model_1preE_newadm_all","ARB_model_1preE_newadm_m","ARB_model_1preE_newadm_f")  # Add last interventions with high-risk groups and pre-emptive isolation
results_epi_ac <- matrix(nrow = 10, ncol = length(model_list))  # Adjust the number of columns based on the number of models/ rows for number of parameters extracted
results_epi_MRSAprev <- matrix(nrow = 366, ncol = length(model_list))
results_epi_MRSAinfe <- matrix(nrow = 366, ncol = length(model_list))
results_epi_MRSAdead <- matrix(nrow = 366, ncol = length(model_list))
results_epi_MRSAinfe_all <- matrix(nrow = 366, ncol = length(model_list))
results_epi_MRSAdead_all <- matrix(nrow = 366, ncol = length(model_list))
results_epi_uncolonised <- matrix(nrow = 366, ncol = length(model_list))

results_econ <- matrix(nrow = 3, ncol = length(model_list))  # Adjust the number of columns based on the number of models/ rows for number of parameters extracted
for (i in 1:length(model_list)) {
  #Run initial conditions first:
  times <- seq(from=0, to=365, by = 1)  # Simulate over a year
  #Staphylococcus aureus params [MRSA/MSSA]####### 
  # ----------------------------------------#
  zeta3_m_p1 <- 0.231  # mortality rate from IMS  , male
  zeta3_f_p1 = zeta3_m_p1*2.07  # mortality rate from IMS  , female
  zeta1_m_p1 = zeta3_m_p1*1.01  # mortality rate from IMR, male
  zeta1_f_p1 = zeta3_m_p1*1.22  # mortality rate from IMR, female
  zeta2_m_p1 = zeta3_m_p1*2.32  # mortality rate from ISR, male
  zeta2_f_p1 = zeta3_m_p1*2.25  # mortality rate from ISR, female
  zeta4_m_p1 = zeta3_m_p1*1.10  # mortality rate from ISS , male
  zeta4_f_p1 = zeta3_m_p1*2.34  # mortality rate from ISS , female
  #(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  kappa_p1 = (35.9-6.2)/(49.5-4.85)
  parameters1 <- c(
    delta1_p1 = 0.0016,  
    delta2_p1 = 0.0016,  
    Disch_U_f_p1 = 1/6,
    Disch_U_m_p1 = 1/6,
    Disch_CR_f_p1 = 1/6,
    Disch_CR_m_p1 = 1/6,
    Disch_CS_f_p1 = 1/6,
    Disch_CS_m_p1 = 1/6,
    mu0_p1 = 0.52,  # % of women among U
    mu1_p1 = 0.52,  # % of women among CR
    mu2_p1 = 0.52,  # % of women among CS
    mu3_p1 = 0.3165,  # % of women among IMR
    mu4_p1 = 0.3889,  # % of women among ISR
    mu5_p1 = 0.3938,  # % of women among IMS
    mu6_p1 = 0.4487,  # % of women among IMS
    #Exposure to anbiotics; treatment amongst susceptible populations
    psi_m_p1 = 0.1474, # % of  individuals exposed  to vancomycin/penicillin among males
    psi_w_p1 = 0.184,  # % of  individuals exposed to vancomycin/penicillin among females
    #Percentage of people under treatment for MRSA decolonisation
    psi_mtr_p1=0.109, #IEAT is 1.35 times higher among CR; hence if treatment is psi_m_p1; psi_mtr_p1=psi_m_p1/1.35
    psi_wtr_p1=0.136,
    #Fitness cost. c reduces the transmission rate among resistant strains [1/unit time] [%].
    c_p1=0.09,
    #Progression to the development of infection from colonisation among CR and CS states. 
    beta1_m_p1 = (1/21)*0.26,  # inverse of LOS plus progression from colonisation to infection among males CR
    beta2_m_p1 = (1/11)*0.099,  # inverse of LOS plus progression from colonisation to infection among males CS
    beta1_f_p1 = (1/29)*0.26, # inverse of LOS plus progression from colonisation to infection among females CR
    beta2_f_p1 = (1/14)*0.099,  # inverse of LOS plus progression from colonisation to infection among females CS
    #Natural clearance of mild and severe infections among CR and CS states, respectively [1/unit time] [%].
    gamma1_p1 = 0.001,  # Natural clearance among mild infections R
    gamma2_p1 = 0.001,  # Natural clearance among severe infections R
    gamma3_p1 = 0.001,  # Natural clearance among mild infections S
    gamma4_p1 = 0.001, # Natural clearance among severe infections S
    #Mean time of infection considering length of hospital stays [1/length of hospital stay]. 
    omega1_d_m_p1 = (1/20), #IMR patients who died, male
    omega1_r_m_p1 = (1/23), #IMR patients who recovered, male
    omega1_d_f_p1 = (1/13), #IMR patients who died, female
    omega1_r_f_p1 = (1/26), #IMR patients who recovered, female
    omega2_d_m_p1 = (1/11),  #ISR patients who died, male
    omega2_r_m_p1 = (1/18),  #ISR patients who recovered, male
    omega2_d_f_p1 = (1/14), #ISR patients who died, female 
    omega2_r_f_p1 = (1/19), #ISR patients who recovered, female
    omega3_d_m_p1 = (1/12), #IMS patients who died, male
    omega3_r_m_p1 = (1/11), #IMS patients who recovered, male
    omega3_d_f_p1 = (1/20), #IMS patients who died, female
    omega3_r_f_p1 = (1/16), #IMS patients who recovered, female
    omega4_d_m_p1 = (1/11), #ISS patients who died, male
    omega4_r_m_p1 = (1/19), #ISS patients who recovered, male
    omega4_d_f_p1 = (1/14), #ISS patients who died, female
    omega4_r_f_p1 = (1/17), #ISS patients who recovered, female
    #Percentage of inpatients with CR or CS, respectively, progressing to severe infection in intensive care units [1/unit time] [%].
    alpha1_m_p1 = 0.4196, #% patients with CR progressing to severe infection, males
    alpha2_m_p1 = 0.3438, #% patients with CS progressing to severe infection, males
    alpha1_f_p1 = 0.3517, #% patients with CR progressing to severe infection, females
    alpha2_f_p1 = 0.2826, #% patients with CS progressing to severe infection, females
    #Progression from mild to severe infection from IMR and IMS, respectively [1/unit time] [%].
    epsilon1_p1 = 0.01, #progression from IMR to ISR
    epsilon2_p1 = 0.01, #progression from IMS to ISS
    #Mortality rates from infection. ζ1 and ζ2 are mortality rates from mild and severe resistant infections, respectively. ζ3 and ζ4 are from mild and severe susceptible infections, respectively [1/unit time] [%].
    zeta3_m_p1 = 0.231,  # mortality rate from IMS  , male
    zeta3_f_p1 = 0.231*2.07,  # mortality rate from IMS  , female
    zeta1_m_p1 = 0.231*1.01,  # mortality rate from IMR, male
    zeta1_f_p1 = 0.231*1.22,  # mortality rate from IMR, female
    zeta2_m_p1 = 0.231*2.32,  # mortality rate from ISR, male
    zeta2_f_p1 = 0.231*2.25,  # mortality rate from ISR, female
    zeta4_m_p1 = 0.231*1.10,  # mortality rate from ISS , male
    zeta4_f_p1 = 0.231*2.34,  # mortality rate from ISS , female
    #Recovery rates from infection, including IMR, ISR, IMS and ISS due to treatment received [1/unit time] [%].
    nu1_m_p1 = (1-zeta1_m_p1), #recovery rates from IMR, males
    nu1_f_p1 = (1-zeta1_f_p1), #recovery rates from IMR, females
    nu2_m_p1 = (1-zeta2_m_p1), #recovery rates from ISR, males
    nu2_f_p1 = (1-zeta2_f_p1), #recovery rates from ISR, females
    nu3_m_p1 = (1-zeta3_m_p1), #recovery rates from IMS, males
    nu3_f_p1 = (1-zeta3_f_p1), #recovery rates from IMS, females
    nu4_m_p1 =  (1-zeta4_m_p1), #recovery rates from ISS, males
    nu4_f_p1 =  (1-zeta4_f_p1),#recovery rates from ISS, females
    #Constant background rate that captures transmission from non-human sources, horizontal transmission, or de novo emergence [1/unit time] [number].
    b_p1 = 0.01,
    # Percentage of people with resistant infections receiving inappropriate empirical antibiotic treatment [1/unit time] [%].
    phi_m_p1 = 0.459,  # 
    phi_f_p1 = 0.413,  # 
    #Factor of burden associated with inappropriate empirical antibiotic treatment and increased ICU admission among resistant infections [1/unit time] [%].
    pi_p1= 1.35,
    #Transmission parameter {update this correspondingly after calibrating it with real data}
    #tau_p1= 0.03461113, #estimated 
    tau_p1=0.2229124,
    #community-acquired infection upon hospital admission rate
    caIha_p1= 0.001,
    ###
    
    #percentage of people tested
    test_p1=0.20, 
    HR_perc1=0.2,
    or_HR_scenar1_a=1.04,
    or_HR_scenarMen_a=1.37,
    ### ### ### ### ###
    #sensitivity chrom_1
    sens_chrom_a=0.826,
    #sensitivity chrom_1
    sens_chrom2_a=0.622, 
    #sensitivity chrom_1
    sens_pcr_a=0.881,
    #turnaround chrom_1
    turn_chrom_a=3,
    #turnaround chrom_1
    turn_chrom2_a=2,
    #turnaround pcr_1
    turn_pcr_a=1,
    #isolation contact precaution transmission reduction
    reduc_conpre_a=0.365,
    #efficiency decolonisation
    eff_decol_a= 0.53,
    #effect on self-infection decolonisation
    eff_decol_selfi_a=0.32,
    #Turnaround decolonisation program in days
    turnaround_decol_a=5,
    ##
    #costs wards
    c_general_ward= 50,
    c_intermediate_ward=92,
    c_icu_ward=218,
    c_decol_1pd=6.5,
    c_isolation=42.3,
    c_chrom=10.2,
    c_chrom2=13.6, 
    c_pcr=33,
    c_bc=16.9,
    #utilities
    u_healthy=0.92,
    u_icu=0.58,
    u_gw=0.64,
    u_recovICU=0.74
  )
  
  
  
  
  # ----------------------------------------#
  #Staphylococcus aureus states [MRSA/MSSA]####### 
  # ----------------------------------------#
  N <- 1000  # Total population size
  # Initial conditions (population sizes in each group)
  U_m10 <- 0.7 * N *(1-0.52)
  CR_m10 <- 0.1194 * N *(1-0.52)
  CS_m10 <-  0.18156* N *(1-0.52)
  IMR_m10 <- 0.26 * CR_m10* (1-0.4196)
  ISR_m10 <- 0.26 * CR_m10* 0.4196
  IMS_m10 <- 0.099 * CS_m10*(1-0.3438)
  ISS_m10 <- 0.099 * CS_m10*0.3438
  RR_m10 <-0
  RS_m10  <-0
  DR_m10  <-0
  DS_m10 <-0
  N_to0 <- 1050
  utility_to0<-0
  cost_to0<-0
  new_admin0<-0
  discharge<-0
  
  U_f10 <- 0.7 * N *0.52
  CR_f10 <- 0.1194 * N *0.52
  CS_f10 <-  0.18156* N *0.52
  IMR_f10 <- 0.26 * CR_f10* (1-0.3517)
  ISR_f10 <- 0.26 * CR_f10*0.3517
  IMS_f10 <- 0.099 * CS_f10*(1-0.2826)
  ISS_f10 <- 0.099 * CS_f10*0.2826
  RR_f10 <-0
  RS_f10  <-0
  DR_f10  <-0
  DS_f10 <-0
  
  N_0m10<-  U_m10 + CR_m10 + CS_m10 + IMR_m10 + ISR_m10 + IMS_m10 + ISS_m10 +  RR_m10 + RS_m10 + DR_m10 +DS_m10
  N_0f10<-  U_f10 + CR_f10 + CS_f10 + IMR_f10 + ISR_f10 + IMS_f10 + ISS_f10 +  RR_f10 + RS_f10 + DR_f10 +DS_f10
  
  state1 <- c(U_m1 = U_m10, CR_m1=CR_m10, CS_m1= CS_m10, IMR_m1= IMR_m10, ISR_m1=ISR_m10, IMS_m1= IMS_m10, ISS_m1= ISS_m10, RR_m1= RR_m10, RS_m1=RS_m10, DR_m1= DR_m10, DS_m1=DS_m10,
              U_f1 = U_f10, CR_f1=CR_f10, CS_f1= CS_f10, IMR_f1= IMR_f10, ISR_f1=ISR_f10, IMS_f1= IMS_f10, ISS_f1= ISS_f10, RR_f1= RR_f10, RS_f1=RS_f10, DR_f1= DR_f10, DS_f1=DS_f10, N_to=N_to0, utility=utility_to0, cost=cost_to0, new_admin=new_admin0, discharge=discharge0)
  
  N_orig1<-N_0m10 + N_0f10
  N_tdif <- N_orig1
  #####
  # Construct the function name from the model_list string
  model_func_name <- get(model_list[i])
  # Solve the ODE using the 'deSolve' package's ode function
  O_solution <- ode(y = state1, times = times, func = model_func_name, parms = parameters1, method = "rk4")
  O_solution2 <- as.data.frame(O_solution)
  O_solution <- as.data.frame(O_solution)
  # Process the ODE solution: sum all values reported in the first 350 rows for your variable of interest
  results_epi_ac[1, i] <- sum(O_solution[1:366, "CR_m1"]) + sum(O_solution[1:366, "CR_f1"]) + sum(O_solution[1:366, "IMR_m1"]) + sum(O_solution[1:366, "IMR_f1"])+ sum(O_solution[1:366, "ISR_m1"])+sum(O_solution[1:366, "ISR_f1"])
  results_epi_ac[2, i] <- sum(O_solution[1:366, "IMR_m1"]) + sum(O_solution[1:366, "IMR_f1"])+ sum(O_solution[1:366, "ISR_m1"])+sum(O_solution[1:366, "ISR_f1"])
  results_epi_ac[3, i] <- (O_solution[366, "DR_f1"]) + (O_solution[366, "DR_m1"]) 
  results_epi_ac[4, i] <- sum(O_solution[1:366, "IMR_m1"]) + sum(O_solution[1:366, "IMR_f1"])+ sum(O_solution[1:366, "ISR_m1"])+sum(O_solution[1:366, "ISR_f1"]) + sum(O_solution[1:366, "IMS_m1"]) + sum(O_solution[1:366, "IMS_f1"])+ sum(O_solution[1:366, "ISS_m1"])+sum(O_solution[1:366, "ISS_f1"])
  results_epi_ac[5, i] <- (O_solution[366, "DR_f1"]) + (O_solution[366, "DR_m1"]) + (O_solution[366, "DS_f1"]) + (O_solution[366, "DS_m1"]) 
  results_epi_ac[6, i] <- sum(O_solution[366, "new_admin"]) 
  results_epi_ac[7, i] <- sum(O_solution[1:366, "U_m1"]) +sum(O_solution[1:366, "U_f1"]) 
  results_epi_ac[8, i] <- sum(O_solution[1:366, "U_m1"]) +sum(O_solution[1:366, "U_f1"]) +sum(O_solution[1:366, "CR_m1"]) +sum(O_solution[1:366, "CR_f1"]) +sum(O_solution[1:366, "CS_m1"]) +sum(O_solution[1:366, "CS_f1"]) 
  results_epi_ac[9, i] <- O_solution[366, "RR_m1"] + O_solution[366, "RR_f1"]+O_solution[366, "RS_m1"] + O_solution[366, "RS_f1"]
  results_epi_ac[10, i] <- O_solution[366, "discharge"]
  
  results_epi_MRSAprev[, i] <- O_solution2$CR_m1 + O_solution2$CR_f1 + O_solution2$IMR_m1 + O_solution2$IMR_f1 + O_solution2$ISR_m1 + O_solution2$ISR_f1
  results_epi_MRSAinfe[, i]  <- O_solution2$IMR_m1 + O_solution2$IMR_f1 + O_solution2$ISR_m1 + O_solution2$ISR_f1
  results_epi_MRSAdead[, i]  <- O_solution2$DR_m1 + O_solution2$DR_f1
  results_epi_MRSAinfe_all[, i]  <-  O_solution2$IMR_m1 + O_solution2$IMR_f1 + O_solution2$ISR_m1 + O_solution2$ISR_f1 + O_solution2$IMS_m1 + O_solution2$IMS_f1 + O_solution2$ISS_m1 + O_solution2$ISS_f1
  results_epi_MRSAdead_all[, i]  <- O_solution2$DR_m1 + O_solution2$DR_f1 + O_solution2$DS_m1 + O_solution2$DS_f1
  
  #Economics
  #Store econ results per strategy
  results_econ[1, i] <- O_solution[366, "cost"]
  results_econ[2, i] <- (results_epi_ac[8, i]*0.92)+(sum(O_solution[1:366, "ISR_m1"]))*0.58+(sum(O_solution[1:366, "ISR_f1"]))*0.58+((sum(O_solution[1:366, "IMR_m1"])+sum(O_solution[1:366, "IMR_f1"]))*0.64)+((sum(O_solution[1:366, "IMS_m1"])+sum(O_solution[1:366, "IMS_f1"]))*0.64)+(sum(O_solution[1:366, "ISS_m1"]))*0.58+(sum(O_solution[1:366, "ISS_f1"]))*0.58 +(O_solution[366, "RR_m1"] + O_solution[366, "RR_f1"]+O_solution[366, "RS_m1"] + O_solution[366, "RS_f1"])*0.92 +results_epi_ac[10, i]*0.92
  results_econ[3, i] <- results_econ[2, i]/O_solution[366,"new_admin"] #check if usage is appropriate
}
#Compute ICER per strategy
results_icer <- matrix(nrow = 13, ncol = length(model_list))
results_icer[1,1] <- 0
results_icer[2,1] <- results_econ[1, 1]
results_icer[3,1] <- results_econ[2, 1]
results_icer[4,1] <- results_epi_ac[1, 1]
results_icer[5,1] <- results_epi_ac[2, 1]
results_icer[6,1] <- results_epi_ac[3, 1]
results_icer[7,1] <- results_epi_ac[4, 1]
results_icer[8,1] <- results_epi_ac[5, 1]
results_icer[9,1] <- results_epi_ac[9, 1]
results_icer[10,1] <- results_epi_ac[8, 1]
results_icer[11,1] <- results_epi_ac[10, 1]
results_icer[12,1] <- results_epi_ac[10, 1]+  results_epi_ac[9, 1]+results_epi_ac[8, 1]+results_epi_ac[4, 1]
results_icer[13,1] <- results_epi_ac[6, 1]
for (i in 2:length(model_list)) {
  results_icer[1, i] <- 0
  results_icer[2, i] <-  results_econ[1, i]
  results_icer[3, i] <-  results_econ[2, i]
  results_icer[4, i]<- results_epi_ac[1, i]
  results_icer[5, i]<- results_epi_ac[2, i]
  results_icer[6, i]<- results_epi_ac[3, i]
  results_icer[7, i]<- results_epi_ac[4, i]
  results_icer[8, i]<- results_epi_ac[5, i]
  results_icer[9, i]<- results_epi_ac[9, i]
  results_icer[10, i]<- results_epi_ac[8, i]
  results_icer[11, i]<- results_epi_ac[10, i]
  results_icer[12, i]<- results_epi_ac[10, i]+  results_epi_ac[9, i]+results_epi_ac[8, i]+results_epi_ac[4, i]
  results_icer[3, i] <- results_econ[2, i] + ifelse((results_icer[12, 1] - results_icer[12, i]) > 0, (results_icer[12, 1] - results_icer[12, i]) * 0.92, 0)
  results_icer[1, i] <- (results_icer[2, i]-results_icer[2, 1])/(results_icer[3, i]-results_icer[3, 1])
  results_icer[13,i] <- results_epi_ac[6, i]
  
}
transposed_icer <- t(results_icer)
transposed_icerdf <- as.data.frame(transposed_icer)
colnames(transposed_icerdf) <- c("ICER", "Costs", "QALYs","ARB colonisation","ARB infections", "ARB deaths", "Total infections","Total deaths","Total Recovered","Total U,CR,CS", "Total discharge","total population")  # Add more names as needed
rownames(transposed_icerdf) <- c("Do-nothing", "Strategy2", "Strategy3", "Strategy4", "Strategy5", "Strategy6", "Strategy7", "Strategy8", "Strategy9", "Strategy10","Strategy11", "Strategy12","Strategy13","Strategy14","Strategy15", "Strategy16")
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/AMR_BSI_BurdenChile/0_Article_AMR Transmission dynamics Chile & LMICs/0_BSIModelling/0_analysis/0_Figures_model")
write.csv(transposed_icerdf, "matrix_resMRSA.csv", row.names = TRUE)

results_epi_acM<- results_epi_ac
results_econM<- results_econ
results_icerM<- results_icer

colonMRSA_a <- matrix(nrow =1, ncol =  length(model_list))
infecMRSA_a <- matrix(nrow =1, ncol =  length(model_list))
DeadMRSA_a <- matrix(nrow =1, ncol =  length(model_list))
Dead_MRSA_tot <- matrix(nrow =1, ncol =  length(model_list))
for (i in 1:length(model_list)) {
  colonMRSA_a[1, i] <- ((results_epi_acM[1, i]-results_epi_acM[2, i])/results_epi_acM[6, i])*100
  infecMRSA_a[1, i] <- (results_epi_acM[2, i]/results_epi_acM[6, i])*100
  DeadMRSA_a[1, i] <-  (results_econM[3, i]/  results_epi_acM[6, i])*100
  Dead_MRSA_tot[1, i] <-  (results_epi_ac[3, i])
}
# Convert matrices to vectors
colonMRSA_a <- as.vector(t(colonMRSA_a))  # Transpose and then convert to vector
infecMRSA_a <- as.vector(t(infecMRSA_a))  # Assuming similar structure
DeadMRSA_a<- as.vector(t(DeadMRSA_a))    # Assuming similar structure
Dead_MRSA_tot<- as.vector(t(Dead_MRSA_tot))



##################################################################################################
##################################################################################################
######CRE#######
##################################################################################################
##################################################################################################
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- ---  MODEL PARAMETERS BELOW --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
#Enterobacterales[CRE/CSE]####### 
# ----------------------------------------#

zeta3_m_p2 = 0.228  # mortality rate from IMS  , male
zeta3_f_p2 = zeta3_m_p2*0.81  # mortality rate from IMS  , female
zeta1_m_p2 = zeta3_m_p2*1.80  # mortality rate from IMR, male
zeta1_f_p2 = zeta3_m_p2*0.55  # mortality rate from IMR, female
zeta2_m_p2 = zeta3_m_p2*1.30  # mortality rate from ISR, male
zeta2_f_p2 = zeta3_m_p2*2.40  # mortality rate from ISR, female
zeta4_m_p2 = zeta3_m_p2*1.62  # mortality rate from ISS , male
zeta4_f_p2 = zeta3_m_p2*2.23  # mortality rate from ISS , female
kappa_p2 = (35.9-6.2)/(49.5-4.85)
parameters2 <- c(
  #clearance -natural- parameters
  delta1_p2 = 0.001,  # Clearance value from CR state
  delta2_p2 = 0.001,  # Clearance value from CS state
  #Dicharge rates from uncolonised and colonised
  Disch_U_f_p2 = 1/6,
  Disch_U_m_p2 = 1/6,
  Disch_CR_f_p2 = 1/6,
  Disch_CR_m_p2 = 1/6,
  Disch_CS_f_p2 = 1/6,
  Disch_CS_m_p2 = 1/6,
  #Proportion of women among specific populations (i.e., U, CR, CS, IMR, ISR, IMS, and ISS) [1/unit time] [%]
  mu0_p2 = 0.52,  # % of women among U
  mu1_p2 = 0.52,  # % of women among CR
  mu2_p2 = 0.52,  # % of women among CS
  mu3_p2 = 0.3394,  # % of women among IMR
  mu4_p2 = 0.3280,  # % of women among ISR
  mu5_p2 = 0.50,  # % of women among IMS
  mu6_p2 = 0.44,  # % of women among IMS
  
  #Exposure to anbiotics 
  psi_m_p2 = 0.2225, # % of  individuals exposed  to vancomycin/penicillin among males
  psi_w_p2 = 0.2026,  # % of  individuals exposed to vancomycin/penicillin among males
  
  #Percentage of people under treatment for CRE decolonisation
  psi_mtr_p2=0.2225/1.02,
  psi_wtr_p2=0.2026/1.02,
  
  #Fitness cost. c reduces the transmission rate among resistant strains [1/unit time] [%].
  c_p2=(1-0.927),
  
  #Progression to the development of infection from colonisation among CR and CS states. 
  beta1_m_p2 = (1/22)*0.213, # inverse of LOS plus progression from colonisation to infection among males CR
  beta2_m_p2 = (1/20)*0.034,  # inverse of LOS plus progression from colonisation to infection among males CS
  beta1_f_p2 = (1/27)*0.213, # inverse of LOS plus progression from colonisation to infection among females CR
  beta2_f_p2 = (1/17)*0.034,  # inverse of LOS plus progression from colonisation to infection among females CS
  
  #Natural clearance of mild and severe infections among CR and CS states, respectively [1/unit time] [%].
  gamma1_p2 = 0.001,  # Natural clearance among mild infections R
  gamma2_p2 = 0.001,  # Natural clearance among severe infections R
  gamma3_p2 = 0.001,  # Natural clearance among mild infections S
  gamma4_p2 = 0.001, # Natural clearance among severe infections S
  
  #Mean time of infection considering length of hospital stays [1/length of hospital stay]. 
  omega1_d_m_p2 =(1/21), #IMR patients who died, male
  omega1_r_m_p2 =(1/26), #IMR patients who recovered, male
  omega1_d_f_p2 =(1/31), #IMR patients who died, female
  omega1_r_f_p2 =(1/30), #IMR patients who recovered, female
  omega2_d_m_p2 =(1/7), #ISR patients who died, male
  omega2_r_m_p2 =(1/20),  #ISR patients who recovered, male
  omega2_d_f_p2 = (1/20), #ISR patients who died, female 
  omega2_r_f_p2 =(1/23), #ISR patients who recovered, female
  omega3_d_m_p2 =(1/12), #IMS patients who died, male
  omega3_r_m_p2 =(1/20), #IMS patients who recovered, male
  omega3_d_f_p2 =(1/10), #IMS patients who died, female
  omega3_r_f_p2 =(1/18), #IMS patients who recovered, female
  omega4_d_m_p2 =(1/11), #ISS patients who died, male
  omega4_r_m_p2 =(1/14), #ISS patients who recovered, male
  omega4_d_f_p2 =(1/9), #ISS patients who died, female
  omega4_r_f_p2 =(1/15), #ISS patients who recovered, female
  
  #Percentage of inpatients with CR or CS, respectively, progressing to severe infection in intensive care units [1/unit time] [%].
  alpha1_f_p2 = 0.4283, #% patients with CR progressing to severe infection, males
  alpha2_f_p2 = 0.3548, #% patients with CS progressing to severe infection, males
  alpha1_m_p2 = 0.4585, #% patients with CR progressing to severe infection, females
  alpha2_m_p2 = 0.3832, #% patients with CS progressing to severe infection, females
  
  #Progression from mild to severe infection from IMR and IMS, respectively [1/unit time] [%].
  epsilon1_p2 = 0.01, #progression from IMR to ISR
  epsilon2_p2 = 0.01, #progression from IMS to ISS
  
  #Mortality rates from infection. ζ1 and ζ2 are mortality rates from mild and severe resistant infections, respectively. ζ3 and ζ4 are from mild and severe susceptible infections, respectively [1/unit time] [%].
  zeta3_m_p2 = 0.228,  # mortality rate from IMS  , male
  zeta3_f_p2 = zeta3_m_p2*0.81,  # mortality rate from IMS  , female
  zeta1_m_p2 = zeta3_m_p2*1.80,  # mortality rate from IMR, male
  zeta1_f_p2 = zeta3_m_p2*0.55,  # mortality rate from IMR, female
  zeta2_m_p2 = zeta3_m_p2*1.30,  # mortality rate from ISR, male
  zeta2_f_p2 = zeta3_m_p2*2.40,  # mortality rate from ISR, female
  zeta4_m_p2 = zeta3_m_p2*1.62,  # mortality rate from ISS , male
  zeta4_f_p2 = zeta3_m_p2*2.23,  # mortality rate from ISS , female
  
  #Recovery rates from infection, including IMR, ISR, IMS and ISS due to treatment received [1/unit time] [%].
  nu1_m_p2 = (1-zeta1_m_p2), #recovery rates from IMR, males
  nu1_f_p2 = (1-zeta1_f_p2), #recovery rates from IMR, females
  nu2_m_p2 = (1-zeta2_m_p2), #recovery rates from ISR, males
  nu2_f_p2 = (1-zeta2_f_p2), #recovery rates from ISR, females
  nu3_m_p2 = (1-zeta3_m_p2), #recovery rates from IMS, males
  nu3_f_p2 = (1-zeta3_f_p2), #recovery rates from IMS, females
  nu4_m_p2 =  (1-zeta4_m_p2), #recovery rates from ISS, males
  nu4_f_p2 =  (1-zeta4_f_p2),#recovery rates from ISS, females
  
  #Constant background rate that captures transmission from non-human sources, horizontal transmission, or de novo emergence [1/unit time] [number].
  b_p2 = 0.01,
  
  # Percentage of people with resistant infections receiving inappropriate empirical antibiotic treatment [1/unit time] [%].
  phi_m_p2 = 0.3782,  # Placeholder value, adjust as needed
  phi_f_p2 = 0.3846,  # Placeholder value, adjust as needed
  
  #Factor of burden associated with inappropriate empirical antibiotic treatment and increased ICU admission among resistant infections [1/unit time] [%].
  pi_p2= 1.02,
  
  #Transmission parameter {update this correspondingly after calibrating it with real data}
  tau_p2= 0.3986551,
  
  #community-acquired infection upon hospital admission rate
  caIha_p2=0.007,
  
  #percentage of people tested
  test_p2=0.20, 
  HR_perc2=0.2,
  or_HR_scenar1_a=1.04,
  or_HR_scenarMen_a=2.27, 
  ### ### ### ### ###
  #sensitivity chrom_1
  sens_chrom_a=0.826,
  #sensitivity chrom_1
  sens_chrom2_a=0.90, 
  #sensitivity chrom_1
  sens_pcr_a=1,
  #turnaround chrom_1
  turn_chrom_a=3,
  #turnaround chrom_1
  turn_chrom2_a=2,
  #turnaround pcr_1
  turn_pcr_a=1,
  #isolation contact precaution transmission reduction
  reduc_conpre_a=0.35,
  #efficiency decolonisation
  eff_decol_a= 0.26, #0.146
  #effect on self-infection decolonisation
  eff_decol_selfi_a=0.041,
  #Turnaround decolonisation program in days
  turnaround_decol_a=7,
  ##
  #costs wards
  c_general_ward= 50,
  c_intermediate_ward=92,
  c_icu_ward=218,
  c_decol_1pd=72.88,
  c_isolation=42.3,
  c_chrom=10.2,
  c_chrom2=13.6, 
  c_pcr=33,
  c_bc=16.9,
  #utilities
  u_healthy=0.92,
  u_icu=0.92-0.34,
  u_gw=0.64,
  u_recovICU=0.74
  
)


######
times <- seq(from=0, to=365, by = 1)  # Simulate over a year
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- ---BASELINE CONDITIONS BELOW--- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
#Enterobacterales[CRE/CSE]####### 
# ----------------------------------#
N <- 1000  # Total population size
# Initial conditions (population sizes in each group)
U_m20 <- 0.44 * N *(1-0.52)
CR_m20 <- 0.1445 * N *(1-0.52)
CS_m20 <-  0.4155* N *(1-0.52)
IMR_m20 <- 0.09 * CR_m20* (1-0.4283)
ISR_m20 <- 0.09 * CR_m20* 0.4283
IMS_m20 <- 0.04 * CS_m20*(1-0.3548)
ISS_m20 <- 0.04 * CS_m20*0.3548
RR_m20 <-0
RS_m20  <-0
DR_m20  <-0
DS_m20 <-0
N_to2<-1050
utility_to0<-0
cost_to0<-0
new_admin0<-0
discharge0<- 0

U_f20 <- 0.44 * N *0.52
CR_f20 <- 0.1445 * N *0.52
CS_f20 <-  0.4155* N *0.52
IMR_f20 <- 0.09 * CR_f20* (1-0.4538)
ISR_f20 <- 0.09 * CR_f20*0.4538
IMS_f20 <- 0.04 * CS_f20*(1-0.3832)
ISS_f20 <- 0.04 * CS_f20*0.3832
RR_f20 <-0
RS_f20  <-0
DR_f20  <-0
DS_f20 <-0

N_0m20<-  U_m20 + CR_m20 + CS_m20 + IMR_m20 + ISR_m20 + IMS_m20 + ISS_m20 +  RR_m20 + RS_m20 + DR_m20 +DS_m20
N_0f20<-  U_f20 + CR_f20 + CS_f20 + IMR_f20 + ISR_f20 + IMS_f20 + ISS_f20 +  RR_f20 + RS_f20 + DR_f20 +DS_f20
N_to0<- N_0m20 + N_0f20

state2 <- c(U_m2 = U_m20, CR_m2=CR_m20, CS_m2= CS_m20, IMR_m2= IMR_m20, ISR_m2=ISR_m20, IMS_m2= IMS_m20, ISS_m2= ISS_m20, RR_m2= RR_m20, RS_m2=RS_m20, DR_m2= DR_m20, DS_m2=DS_m20,
            U_f2 = U_f20, CR_f2=CR_f20, CS_f2= CS_f20, IMR_f2= IMR_f20, ISR_f2=ISR_f20, IMS_f2= IMS_f20, ISS_f2= ISS_f20, RR_f2= RR_f20, RS_f2=RS_f20, DR_f2= DR_f20, DS_f2=DS_f20, N_to2=N_to0, utility=utility_to0, cost=cost_to0, new_admin=new_admin0, discharge=discharge0)
N_orig2<-N_0m20 + N_0f20
N_tdif <- N_orig2
Nt2_spec2<-1030
#####
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- MODEL & EQUATIONS BELOW --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# Define a function for the differential equations
#N original baseline conditions
ARB_model_2ch_do_nothing <- function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_1 <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_1 <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_2 <- (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2)+ b_p2*(r_v22)) + ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)
  
  #Interventions
  interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a+1))
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+(INF_U_m_p2)-(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+INF_CR_m_p2+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2)
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2)
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+INF_CR_f_p2+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2)
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2)
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
#I.1. test+ treatment decolonisation, all new admissions 
ARB_model_2ch_td_newadm <- function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_2 <- (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  FOC_cr_2+FOC_cs_2
  
  #Interventions
  interv_inf_Rpd1<- (sens_chrom_a*eff_decol_a)*(1/(turn_chrom_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a)/(turn_chrom_a)
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+(INF_U_m_p2 + INF_CR_m_p2*(interv_inf_Rpd1)+ INF_IMR_m_p2*(interv_inf_Rpd1) + INF_ISR_m_p2*(interv_inf_Rpd1))-(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+(INF_CR_m_p2*(1-interv_inf_Rpd1) + interv_inf_reductPr*INF_IMR_m_p2 +interv_inf_reductPr*INF_ISR_m_p2)+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+ INF_U_f_p2+ ((INF_CR_f_p2*(interv_inf_Rpd1)) +(INF_IMR_f_p2*(interv_inf_Rpd1)) + (INF_ISR_f_p2*(interv_inf_Rpd1)))-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+(INF_CR_f_p2*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_f_p2 +interv_inf_reductPr*INF_ISR_f_p2)+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (influx_nonARB)*(c_chrom) + (influx_ARB)*(c_chrom+c_decol_1pd) +c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
ARB_model_2ch2_td_newadm <-function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<- state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_2 <- (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2)+ b_p2*(r_v22)) + ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)
  
  #Interventions
  interv_inf_Rpd1<- (sens_chrom2_a*eff_decol_a)*(1/(turn_chrom2_a+turnaround_decol_a))
  interv_inf_reductPr <-((eff_decol_selfi_a)/(turn_chrom2_a))
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+INF_U_m_p2 + ((INF_CR_m_p2*(interv_inf_Rpd1)) +(INF_IMR_m_p2*(interv_inf_Rpd1)) + (INF_ISR_m_p2*(interv_inf_Rpd1)))-(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+(INF_CR_m_p2*(1-interv_inf_Rpd1)) +(interv_inf_reductPr*INF_IMR_m_p2 +interv_inf_reductPr*INF_ISR_m_p2)+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+ INF_U_f_p2+((INF_CR_f_p2*(interv_inf_Rpd1)) +(INF_IMR_f_p2*(interv_inf_Rpd1)) + (INF_ISR_f_p2*(interv_inf_Rpd1)))-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+(INF_CR_f_p2*(1-interv_inf_Rpd1)) +(interv_inf_reductPr*INF_IMR_f_p2 +interv_inf_reductPr*INF_ISR_f_p2)+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (influx_nonARB)*(c_chrom2) + (influx_ARB)*(c_chrom2+c_decol_1pd) +c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
ARB_model_2pcr_td_newadm <-function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_1 <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_1 <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_2 <- (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  FOC_cr_2+FOC_cs_2
  
  #Interventions
  interv_inf_Rpd1<- (sens_pcr_a*eff_decol_a)*(1/(turn_pcr_a+turnaround_decol_a))
  interv_inf_reductPr <-((eff_decol_selfi_a)/(turn_pcr_a))
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+(INF_U_m_p2) + ((INF_CR_m_p2*(interv_inf_Rpd1)) +(INF_IMR_m_p2*(interv_inf_Rpd1)) + (INF_ISR_m_p2*(interv_inf_Rpd1)))-(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+(INF_CR_m_p2*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_m_p2 +interv_inf_reductPr*INF_ISR_m_p2)+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2) + ((INF_CR_f_p2*(interv_inf_Rpd1)) +(INF_IMR_f_p2*(interv_inf_Rpd1)) + (INF_ISR_f_p2*(interv_inf_Rpd1)))-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+(INF_CR_f_p2*(1-interv_inf_Rpd1)) +(interv_inf_reductPr*INF_IMR_f_p2 +interv_inf_reductPr*INF_ISR_f_p2)+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2*(1-interv_inf_Rpd1-interv_inf_reductPr)) 
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (influx_nonARB)*(c_pcr) + (influx_ARB)*(c_pcr+c_decol_1pd) +c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
#I.2. test+treatment isolation, all new admissions 
ARB_model_2ch_tiso_newadm <-function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge <- state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_1 <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_1 <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  #Interventions
  interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a))
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_2 <- (((tau_p2*(1-c_p2)*(1-(sens_chrom_a*(1/(turn_chrom_a))*reduc_conpre_a))*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  FOC_cr_2+FOC_cs_2
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+(INF_U_m_p2)-(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+INF_CR_m_p2+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2)
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2)
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+INF_CR_f_p2+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2)
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2)
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (influx_nonARB)*(c_chrom) + (influx_ARB)*(c_chrom+ c_isolation)+c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
ARB_model_2ch2_tiso_newadm <-function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<- state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_1 <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_1 <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  #Interventions
  interv_inf_Rpd1<- sens_chrom2_a*eff_decol_a*(1/(turn_chrom2_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom2_a))
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_2 <- (((tau_p2*(1-c_p2)*(1-(sens_chrom2_a*(1/(turn_chrom2_a))*reduc_conpre_a))*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  FOC_cr_2+FOC_cs_2
  
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+(INF_U_m_p2)-(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+INF_CR_m_p2+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2)
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2)
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+INF_CR_f_p2+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2)
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2)
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (influx_nonARB)*(c_chrom2) + (influx_ARB)*(c_chrom2+ c_isolation)+c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
ARB_model_2pcr_tiso_newadm <-function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_1 <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_1 <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  
  #Interventions
  interv_inf_Rpd1<- sens_pcr_a*eff_decol_a*(1/(turn_pcr_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_pcr_a))
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_2 <- (((tau_p2*(1-c_p2)*(1-(sens_pcr_a*(1/(turn_pcr_a))*reduc_conpre_a))*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  FOC_cr_2+FOC_cs_2
  
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+(INF_U_m_p2)-(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+INF_CR_m_p2+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2)
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2)
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+INF_CR_f_p2+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2)
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2)
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (influx_nonARB)*(c_pcr) + (influx_ARB)*(c_pcr+ c_isolation)+c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
#II.1. treatment decolonisation, high-risk patients: males
ARB_model_2ch_td_newadmHR_m <-function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_1 <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_1 <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_2 <- (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2)+ b_p2*(r_v22)) + ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)
  
  #Interventions
  interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a))
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  Influx_men <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2 +INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_r <-  INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_nor <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2
  Influx_women <-INF_U_f_p2 + INF_CS_f_p2+ INF_IMS_f_p2 + INF_ISS_f_p2+INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_r <-  INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_nor <-  INF_U_f_p2  + INF_CS_f_p2 + INF_IMS_f_p2 + INF_ISS_f_p2
  
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+INF_U_m_p2+  (INF_CR_m_p2*(interv_inf_Rpd1)) +(INF_IMR_m_p2*(interv_inf_Rpd1))+ (INF_ISR_m_p2*(interv_inf_Rpd1)) -(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)+(INF_CR_m_p2)*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_m_p2 +interv_inf_reductPr*INF_ISR_m_p2
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2)*(1-interv_inf_Rpd1-interv_inf_reductPr) 
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2)*(1-interv_inf_Rpd1-interv_inf_reductPr) 
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+INF_CR_f_p2+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2)
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2)
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (Influx_men_nor)*(c_chrom) + (Influx_men_r)*(c_chrom+c_decol_1pd) + c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
ARB_model_2ch2_td_newadmHR_m <-function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<- state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_1 <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_1 <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_2 <- (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2)+ b_p2*(r_v22)) + ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)
  
  #Interventions
  interv_inf_Rpd1<- sens_chrom2_a*eff_decol_a*(1/(turn_chrom2_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom2_a))
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  Influx_men <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2 +INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_r <-  INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_nor <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2
  Influx_women <-INF_U_f_p2 + INF_CS_f_p2+ INF_IMS_f_p2 + INF_ISS_f_p2+INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_r <-  INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_nor <-  INF_U_f_p2  + INF_CS_f_p2 + INF_IMS_f_p2 + INF_ISS_f_p2
  
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+INF_U_m_p2+  (INF_CR_m_p2*(interv_inf_Rpd1)) +(INF_IMR_m_p2*(interv_inf_Rpd1))+ (INF_ISR_m_p2*(interv_inf_Rpd1)) -(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)+(INF_CR_m_p2)*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_m_p2 +interv_inf_reductPr*INF_ISR_m_p2
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2)*(1-interv_inf_Rpd1-interv_inf_reductPr) 
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2)*(1-interv_inf_Rpd1-interv_inf_reductPr) 
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+INF_CR_f_p2+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2)
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2)
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (Influx_men_nor)*(c_chrom2) + (Influx_men_r)*(c_chrom2+c_decol_1pd) + c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
ARB_model_2_pcr_td_newadmHR_m <-function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_1 <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_1 <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_2 <- (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2)+ b_p2*(r_v22)) + ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)
  
  #Interventions
  interv_inf_Rpd1<- sens_pcr_a*eff_decol_a*(1/(turn_pcr_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_pcr_a))
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  Influx_men <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2 +INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_r <-  INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_nor <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2
  Influx_women <-INF_U_f_p2 + INF_CS_f_p2+ INF_IMS_f_p2 + INF_ISS_f_p2+INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_r <-  INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_nor <-  INF_U_f_p2  + INF_CS_f_p2 + INF_IMS_f_p2 + INF_ISS_f_p2
  
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+INF_U_m_p2+  (INF_CR_m_p2*(interv_inf_Rpd1)) +(INF_IMR_m_p2*(interv_inf_Rpd1))+ (INF_ISR_m_p2*(interv_inf_Rpd1)) -(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)+(INF_CR_m_p2)*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_m_p2 +interv_inf_reductPr*INF_ISR_m_p2
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2)*(1-interv_inf_Rpd1-interv_inf_reductPr) 
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2)*(1-interv_inf_Rpd1-interv_inf_reductPr) 
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+INF_CR_f_p2+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2)
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2)
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (Influx_men_nor)*(c_pcr) + (Influx_men_r)*(c_pcr+c_decol_1pd) + c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
#II.2. treatment decolonisation, high-risk patients: females  
ARB_model_2ch_td_newadmHR_f <- function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge <-state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_1 <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_1 <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_2 <- (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2)+ b_p2*(r_v22)) + ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)
  
  #Interventions
  interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a))
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  Influx_men <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2 +INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_r <-  INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_nor <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2
  Influx_women <-INF_U_f_p2 + INF_CS_f_p2+ INF_IMS_f_p2 + INF_ISS_f_p2+INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_r <-  INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_nor <-  INF_U_f_p2  + INF_CS_f_p2 + INF_IMS_f_p2 + INF_ISS_f_p2
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+(INF_U_m_p2)-(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+INF_CR_m_p2+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2)
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2)
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)+(INF_IMR_f_p2*(interv_inf_Rpd1)) + (INF_ISR_f_p2*(interv_inf_Rpd1))+INF_CR_f_p2*(interv_inf_Rpd1)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)+(INF_CR_f_p2*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_f_p2 +interv_inf_reductPr*INF_ISR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (Influx_women_nor)*(c_chrom) + (Influx_women_r)*(c_chrom+c_decol_1pd) + c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
ARB_model_2ch2_td_newadmHR_f <- function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<-state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_1 <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_1 <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_2 <- (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2)+ b_p2*(r_v22)) + ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)
  
  #Interventions
  interv_inf_Rpd1<- sens_chrom2_a*eff_decol_a*(1/(turn_chrom2_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom2_a))
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  Influx_men <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2 +INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_r <-  INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_nor <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2
  Influx_women <-INF_U_f_p2 + INF_CS_f_p2+ INF_IMS_f_p2 + INF_ISS_f_p2+INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_r <-  INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_nor <-  INF_U_f_p2  + INF_CS_f_p2 + INF_IMS_f_p2 + INF_ISS_f_p2
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+(INF_U_m_p2)-(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+INF_CR_m_p2+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2)
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2)
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)+(INF_IMR_f_p2*(interv_inf_Rpd1)) + (INF_ISR_f_p2*(interv_inf_Rpd1)) +INF_CR_f_p2*(interv_inf_Rpd1)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)+(INF_CR_f_p2*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_f_p2 +interv_inf_reductPr*INF_ISR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (Influx_women_nor)*(c_chrom2) + (Influx_women_r)*(c_chrom2+c_decol_1pd) + c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
ARB_model_2_pcr_td_newadmHR_f <-function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<- state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_1 <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_1 <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  FOC_cr_2 <- (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  (((tau_p2*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2)+ b_p2*(r_v22)) + ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)
  
  #Interventions
  interv_inf_Rpd1<- sens_pcr_a*eff_decol_a*(1/(turn_pcr_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_pcr_a))
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  Influx_men <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2 +INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_r <-  INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_nor <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2
  Influx_women <-INF_U_f_p2 + INF_CS_f_p2+ INF_IMS_f_p2 + INF_ISS_f_p2+INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_r <-  INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_nor <-  INF_U_f_p2  + INF_CS_f_p2 + INF_IMS_f_p2 + INF_ISS_f_p2
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+(INF_U_m_p2)-(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+INF_CR_m_p2+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2)
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2)
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)+(INF_IMR_f_p2*(interv_inf_Rpd1)) + (INF_ISR_f_p2*(interv_inf_Rpd1)) +INF_CR_f_p2*(interv_inf_Rpd1)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)+(INF_CR_f_p2*(1-interv_inf_Rpd1) +interv_inf_reductPr*INF_IMR_f_p2 +interv_inf_reductPr*INF_ISR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2)*(1-interv_inf_Rpd1-interv_inf_reductPr)
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (Influx_women_nor)*(c_pcr) + (Influx_women_r)*(c_pcr+c_decol_1pd) + c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
#III.1 Pre-emptive isolation of all new admissions [NO TEST]
ARB_model_2preE_newadm_all <- function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge <-state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_1 <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_1 <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  #RESCALING TRANSMISSION PARAMETER
  #tau_p1_rs <- (75.14+69.36)*tau_p2 /((1-0.6069)*(75.14)+0.6069*(69.36)) #(CR_m+CR_f1)*tau_p1 /((1-0.6069)*(CR_f1)+0.6069*(CR_m))
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions  0.3*Beta*women*uncolonised/(Nt) + 0.7*(1-Clevel)*Beta*(men)*uncolonised/Nt  Clevel= Coverage*efficacy*OR
  FOC_cr_2 <- (((tau_p2*(1-reduc_conpre_a)*(1-c_p2)*((CR_f2+CR_m2+IMR_m2+IMR_f2+ISR_m2+ISR_f2)*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22))
  FOC_cs_2 <- ((tau_p2*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  FOC_cs_2+FOC_cr_2
  
  #Interventions
  interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a+1))
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+(INF_U_m_p2)-(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+INF_CR_m_p2+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2)
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2)
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+INF_CR_f_p2+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2)
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2)
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (influx_nonARB + influx_ARB)*(c_isolation)+c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
ARB_model_2preE_newadm_m <- function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge <-state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_a <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_a <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  #RESCALING TRANSMISSION PARAMETER
  tau_p2_rs <- (75.14+69.36)*tau_p2 /((1-0.48)*(75.14)+0.48*(69.36)) #(CR_m+CR_f1)*tau_p1 /((1-0.6069)*(CR_f1)+0.6069*(CR_m))
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions  0.3*Beta*women*uncolonised/(Nt) + 0.7*(1-Clevel)*Beta*(men)*uncolonised/Nt  Clevel= Coverage*efficacy*OR
  FOC_cr_2 <- ((((1-0.48)*tau_p2_rs*(1-c_p2)*((CR_f2+IMR_f2+ISR_f2))*(U_f2+U_m2)))/Nt2_spec2)+(((0.48*tau_p2_rs*(1-c_p2)*((CR_m2+IMR_m2+ISR_m2)*(1-reduc_conpre_a*or_HR_scenarMen_a))*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22)
  FOC_cs_2 <- ((tau_p2*((CS_m2+IMS_f2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  + ((tau_p2*((CS_m2+IMS_m2+ISS_m2)*(1-reduc_conpre_a)*(U_f2+U_m2)))/Nt2_spec2)  
  FOC_u_2 <-  FOC_cs_2+FOC_cr_2
  
  #Interventions
  interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a+1))
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  Influx_men <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2 +INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_r <-  INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_nor <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2
  Influx_women <-INF_U_f_p2 + INF_CS_f_p2+ INF_IMS_f_p2 + INF_ISS_f_p2+INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_r <-  INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_nor <-  INF_U_f_p2  + INF_CS_f_p2 + INF_IMS_f_p2 + INF_ISS_f_p2
  
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+(INF_U_m_p2)-(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+INF_CR_m_p2+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2)
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2)
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+INF_CR_f_p2+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2)
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2)
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (Influx_men)*(c_isolation)+c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}
ARB_model_2preE_newadm_f <- function(times, state, parms) {
  ## Define variables
  
  # Men 
  U_m2 <- state["U_m2"]
  CR_m2 <- state["CR_m2"]
  CS_m2 <- state["CS_m2"]
  IMR_m2 <- state["IMR_m2"]
  ISR_m2 <- state["ISR_m2"]
  IMS_m2 <- state["IMS_m2"]
  ISS_m2 <- state["ISS_m2"]
  RR_m2 <- state["RR_m2"]
  RS_m2 <- state["RS_m2"]
  DR_m2 <- state["DR_m2"]
  DS_m2 <- state["DS_m2"]
  
  N1_2 <- U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 + RR_m2 + RS_m2 + DR_m2 + DS_m2
  
  # Women   
  U_f2 <- state["U_f2"]
  CR_f2 <- state["CR_f2"]
  CS_f2 <- state["CS_f2"]
  IMR_f2 <- state["IMR_f2"]
  ISR_f2 <- state["ISR_f2"]
  IMS_f2 <- state["IMS_f2"]
  ISS_f2 <- state["ISS_f2"]
  RR_f2 <- state["RR_f2"]
  RS_f2 <- state["RS_f2"]
  DR_f2 <- state["DR_f2"]
  DS_f2 <- state["DS_f2"]
  N_to2 <- state["N_to"]
  utility <- state["utility"]
  cost <- state["cost"]
  new_admin<-state["new_admin"]
  discharge<- state["discharge"]
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 =  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  
  # # # # # # #
  
  #Extract parameters
  delta1_p2<- parms["delta1_p2"]
  delta2_p2<- parms["delta2_p2"]
  Disch_U_f_p2<-parms["Disch_U_f_p2"]
  Disch_U_m_p2<- parms["Disch_U_m_p2"] 
  Disch_CR_f_p2<- parms["Disch_CR_f_p2"] 
  Disch_CR_m_p2<-parms["Disch_CR_m_p2"] 
  Disch_CS_f_p2<-parms["Disch_CS_f_p2"] 
  Disch_CS_m_p2 <-parms["Disch_CS_m_p2"]
  mu0_p2<- parms["mu0_p2"]
  mu1_p2<- parms["mu1_p2"]
  mu2_p2<- parms["mu2_p2"]
  mu3_p2<- parms["mu3_p2"]
  mu4_p2<- parms["mu4_p2"]
  mu5_p2<- parms["mu5_p2"]
  mu6_p2<- parms["mu6_p2"]
  psi_m_p2<- parms["psi_m_p2"]
  psi_w_p2<- parms["psi_w_p2"]
  c_p2<- parms["c_p2"]
  beta1_m_p2<- parms["beta1_m_p2"]
  beta2_m_p2 <- parms["beta2_m_p2"]
  beta1_f_p2<- parms["beta1_f_p2"]
  beta2_f_p2<- parms["beta2_f_p2"]
  gamma1_p2<- parms["gamma1_p2"]
  gamma2_p2<- parms["gamma2_p2"]
  gamma3_p2<- parms["gamma3_p2"]
  gamma4_p2<- parms["gamma4_p2"]
  omega1_d_m_p2<- parms["omega1_d_m_p2"]
  omega1_r_m_p2<- parms["omega1_r_m_p2"]
  omega1_d_f_p2<- parms["omega1_d_f_p2"]
  omega1_r_f_p2<- parms["omega1_r_f_p2"]
  omega2_d_m_p2<- parms["omega2_d_m_p2"]
  omega2_r_m_p2<- parms["omega2_r_m_p2"]
  omega2_d_f_p2<- parms["omega2_d_f_p2"]
  omega2_r_f_p2<- parms["omega2_r_f_p2"]
  omega3_d_m_p2<- parms["omega3_d_m_p2"]
  omega3_r_m_p2<- parms["omega3_r_m_p2"]
  omega3_d_f_p2<- parms["omega3_d_f_p2"]
  omega3_r_f_p2<- parms["omega3_r_f_p2"]
  omega4_d_m_p2<- parms["omega4_d_m_p2"]
  omega4_r_m_p2<- parms["omega4_r_m_p2"]
  omega4_d_f_p2<- parms["omega4_d_f_p2"]
  omega4_r_f_p2<- parms["omega4_r_f_p2"]
  alpha1_m_p2<- parms["alpha1_m_p2"]
  alpha2_m_p2<- parms["alpha2_m_p2"]
  alpha1_f_p2<- parms["alpha1_f_p2"]
  alpha2_f_p2<- parms["alpha2_f_p2"]
  epsilon1_p2<- parms["epsilon1_p2"]
  epsilon2_p2<- parms["epsilon2_p2"]
  zeta3_m_p2<- parms["zeta3_m_p2"]
  zeta3_f_p2<- parms["zeta3_f_p2"]
  zeta1_m_p2<- parms["zeta1_m_p2"]
  zeta1_f_p2<- parms["zeta1_f_p2"]
  zeta2_m_p2<- parms["zeta2_m_p2"]
  zeta2_f_p2<- parms["zeta2_f_p2"]
  zeta4_m_p2<- parms["zeta4_m_p2"]
  zeta4_f_p2<- parms["zeta4_f_p2"]
  nu1_m_p2<- parms["nu1_m_p2"]
  nu1_f_p2<- parms["nu1_f_p2"]
  nu2_m_p2 <- parms["nu2_m_p2"]
  nu2_f_p2<- parms["nu2_f_p2"]
  nu3_m_p2<- parms["nu3_m_p2"]
  nu3_f_p2<- parms["nu3_f_p2"]
  nu4_m_p2<- parms["nu4_m_p2"]
  nu4_f_p2<- parms["nu4_f_p2"]
  b_p2<- parms["b_p2"]
  phi_m_p2<- parms["phi_m_p2"]
  phi_f_p2<- parms["phi_f_p2"]
  pi_p2<- parms["pi_p2"]
  tau_p2<-parms["tau_p2"]
  caIha_p2<-parms["caIha_p2"]
  psi_mtr_p2 <-parms["psi_mtr_p2"]
  psi_wtr_p2 <-parms["psi_wtr_p2"]
  test_p2 <-parms["test_p2"]
  test_p2<-parms["test_p2"] 
  or_HR_scenar1_1 <- parms["or_HR_scenar1_a"]
  or_HR_scenarMen_1 <- parms["or_HR_scenarMen_a"]
  #sensitivity chrom_1
  sens_chrom_a <- parms["sens_chrom_a"]
  #sensitivity chrom_1
  sens_chrom2_a <- parms["sens_chrom2_a"]
  #sensitivity chrom_1
  sens_pcr_a <- parms["sens_pcr_a"]
  #turnaround chrom_1
  turn_chrom_a <- parms["turn_chrom_a"]  
  #turnaround chrom_1
  turn_chrom2_a <- parms["turn_chrom2_a"]
  #turnaround pcr_1
  turn_pcr_a <- parms["turn_pcr_a"]
  #isolation contact precaution transmission reduction
  reduc_conpre_a <- parms["reduc_conpre_a"]
  #efficiency decolonisation
  eff_decol_a <- parms["eff_decol_a"]
  #effect on self-infection decolonisation
  eff_decol_selfi_a <- parms["eff_decol_selfi_a"]
  #Turnaround decolonisation program in days
  turnaround_decol_a <- parms["turnaround_decol_a"] 
  ##
  #cost hospital wards
  c_general_ward <- parms["c_general_ward"]
  c_intermediate_ward <- parms["c_intermediate_ward"]
  c_icu_ward <- parms["c_icu_ward"]
  c_decol_1pd <- parms["c_decol_1pd"]
  c_isolation <- parms["c_isolation"]
  c_chrom <- parms["c_chrom"]
  c_chrom2 <- parms["c_chrom2"]
  c_pcr <- parms["c_pcr"]
  c_bc <- parms["c_bc"]
  #utilities
  u_healthy <- parms["u_healthy"]
  u_icu <- parms["u_icu"]
  u_gw <- parms["u_gw"]
  u_recovICU <- parms["u_recovICU"]
  
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  
  #Influx of populations
  INF_U_f_p2 <- (1050- Nt2_spec2)*0.44*mu0_p2
  INF_U_m_p2 <- (1050- Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050- Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2<- (1050- Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2<- (1050- Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050- Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050- Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050- Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
  #Prevalence of CRE
  P1_t2 <- (CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2)/(CR_f2 + IMR_f2 + ISR_f2 + CR_m2 + IMR_m2 + ISR_m2 + CS_f2 + IMS_f2 + ISS_f2 + CS_m2 + IMS_m2 + ISS_m2)
  
  #Random value for competing transmissions
  ra_v2 <- runif(1, min = 0.00, max = 0.01)
  alpha12 <- 0.5
  beta12 <- (0.5)
  r_v22 <- rbeta(1, alpha12, beta12)
  h_ieat1_p2 <- (alpha1_m_p2)/((pi_p2*phi_m_p2)+(1-phi_m_p2))
  h_ieat2_p2 <- (alpha1_f_p2)/((pi_p2*phi_f_p2)+(1-phi_f_p2))
  N_to2<- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions 
  #RESCALING TRANSMISSION PARAMETER
  tau_p2_rs <- (75.14+69.36)*tau_p2 /((1-0.48)*(75.14)+0.48*(69.36)) #(CR_m+CR_f1)*tau_p1 /((1-0.6069)*(CR_f1)+0.6069*(CR_m))
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions  0.3*Beta*women*uncolonised/(Nt) + 0.7*(1-Clevel)*Beta*(men)*uncolonised/Nt  Clevel= Coverage*efficacy*OR
  FOC_cr_2 <- ((((1-0.48)*tau_p2_rs*(1-c_p2)*((CR_f2+IMR_f2+ISR_f2))*(1-reduc_conpre_a)*(U_f2+U_m2)))/Nt2_spec2)+(((0.48*tau_p2_rs*(1-c_p2)*((CR_m2+IMR_m2+ISR_m2))*(U_f2+U_m2)))/Nt2_spec2) + b_p2*(r_v22)
  FOC_cs_2 <- ((tau_p2*((CS_m2+IMS_m2+ISS_m2)*(U_f2+U_m2)))/Nt2_spec2)  + ((tau_p2*((CS_f2+IMS_f2+ISS_f2)*(1-reduc_conpre_a)*(U_f2+U_m2)))/Nt2_spec2)
  FOC_u_2 <-  FOC_cs_2+FOC_cr_2
  
  #Interventions
  interv_inf_Rpd1<- sens_chrom_a*eff_decol_a*(1/(turn_chrom_a+turnaround_decol_a))
  interv_inf_reductPr <-(eff_decol_selfi_a/(turn_chrom_a+1))
  
  #INFLUX FOR INTERVENTIONS!
  influx_nonARB<- INF_U_f_p2 + INF_U_m_p2 + INF_CS_f_p2 + INF_CS_m_p2 + INF_IMS_f_p2 + INF_IMS_m_p2 + INF_ISS_f_p2 + INF_ISS_m_p2
  influx_ARB<- INF_CR_f_p2 + INF_CR_m_p2+  INF_IMR_f_p2 + INF_IMR_m_p2 + INF_ISR_f_p2 + INF_ISR_m_p2
  Influx_men <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2 +INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_r <-  INF_CR_m_p2+INF_IMR_m_p2 +INF_ISR_m_p2
  Influx_men_nor <-  INF_U_m_p2  + INF_CS_m_p2 + INF_IMS_m_p2 + INF_ISS_m_p2
  Influx_women <-INF_U_f_p2 + INF_CS_f_p2+ INF_IMS_f_p2 + INF_ISS_f_p2+INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_r <-  INF_CR_f_p2 +INF_IMR_f_p2 +INF_ISR_f_p2
  Influx_women_nor <-  INF_U_f_p2  + INF_CS_f_p2 + INF_IMS_f_p2 + INF_ISS_f_p2
  
  # DEFINITION OF THE DIFFERENTIAL EQUATIONS
  
  dU_m2 <-   (delta1_p2*CR_m2)+(delta2_p2*CS_m2)+(INF_U_m_p2)-(FOC_u_2*(1-mu0_p2))-(U_m2*Disch_U_m_p2)+(psi_m_p2*CS_m2)+(psi_mtr_p2*CR_m2)
  dCR_m2 <- -(delta1_p2*CR_m2)-(beta1_m_p2*CR_m2)-(psi_mtr_p2*CR_m2)+(gamma1_p2*IMR_m2)+(gamma2_p2*ISR_m2)+INF_CR_m_p2+((1-mu1_p2)*(FOC_cr_2))-(CR_m2*Disch_CR_m_p2)
  dCS_m2 <- -(delta2_p2*CS_m2)-(beta2_m_p2*CS_m2)-(psi_m_p2*CS_m2)  +(gamma3_p2*IMS_m2)+(gamma4_p2*ISS_m2)+INF_CS_m_p2+((1-mu2_p2)*(FOC_cs_2))-(CS_m2*Disch_CS_m_p2)
  dIMR_m2 <- ((beta1_m_p2*CR_m2)*(1-alpha1_m_p2))-(gamma1_p2*IMR_m2)-(omega1_r_m_p2*nu1_m_p2*IMR_m2)-(epsilon1_p2*IMR_m2)-(omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(INF_IMR_m_p2)
  dISR_m2 <- (beta1_m_p2*CR_m2*alpha1_m_p2)      -(gamma2_p2*ISR_m2)-(omega2_r_m_p2*nu2_m_p2*ISR_m2)+(epsilon1_p2*IMR_m2)-(omega2_d_m_p2*zeta2_m_p2*ISR_m2)+(INF_ISR_m_p2)
  dIMS_m2 <- (beta2_m_p2*CS_m2*(1-alpha2_m_p2))-(gamma3_p2*IMS_m2)-(omega3_r_m_p2*nu3_m_p2*IMS_m2)-(epsilon2_p2*ISS_m2)-(omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(INF_IMS_m_p2)
  dISS_m2 <- (beta2_m_p2*CS_m2*(alpha2_m_p2))-(gamma4_p2*ISS_m2)-(omega4_r_m_p2*nu4_m_p2*ISS_m2)+(epsilon2_p2*ISS_m2)-(omega4_d_m_p2*zeta4_m_p2*ISS_m2)+(INF_ISS_m_p2)
  dRR_m2 <-  (omega1_r_m_p2*nu1_m_p2*IMR_m2)+(omega2_r_m_p2*nu2_m_p2*ISR_m2)
  dRS_m2 <-  (omega3_r_m_p2*nu3_m_p2*IMS_m2)+(omega4_r_m_p2*nu4_m_p2*ISS_m2)
  dDR_m2 <-  (omega1_d_m_p2*zeta1_m_p2*IMR_m2)+(omega2_d_m_p2*zeta2_m_p2*ISR_m2)
  dDS_m2 <-  (omega3_d_m_p2*zeta3_m_p2*IMS_m2)+(omega4_d_m_p2*zeta4_m_p2*ISS_m2)
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_f2)+(psi_wtr_p2*CR_f2)
  dCR_f2 <- -(delta1_p2*CR_f2)-(beta1_f_p2*CR_f2)-(psi_wtr_p2*CR_f2)+(gamma1_p2*IMR_f2)+(gamma2_p2*ISR_f2)+INF_CR_f_p2+((mu1_p2)*(FOC_cr_2))-(CR_f2*Disch_CR_f_p2)
  dCS_f2<-  -(delta2_p2*CS_f2)-(beta2_f_p2*CS_f2)-(psi_w_p2*CS_f2)  +(gamma3_p2*IMS_f2)+(gamma4_p2*ISS_f2)+INF_CS_f_p2+((mu2_p2)*(FOC_cs_2))-(CS_f2*Disch_CS_f_p2)
  dIMR_f2 <- ((beta1_f_p2*CR_f2)*(1-alpha1_f_p2))-(gamma1_p2*IMR_f2)-(omega1_r_f_p2*nu1_f_p2*IMR_f2)-(epsilon1_p2*IMR_f2)-(omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(INF_IMR_f_p2)
  dISR_f2 <- (beta1_f_p2*CR_f2*alpha1_f_p2)      -(gamma2_p2*ISR_f2)-(omega2_r_f_p2*nu2_f_p2*ISR_f2)+(epsilon1_p2*IMR_f2)-(omega2_d_f_p2*zeta2_f_p2*ISR_f2)+(INF_ISR_f_p2)
  dIMS_f2 <- (beta2_f_p2*CS_f2*(1-alpha2_f_p2))  -(gamma3_p2*IMS_f2)-(omega3_r_f_p2*nu3_f_p2*IMS_f2)-(epsilon2_p2*ISS_f2)-(omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(INF_IMS_f_p2)
  dISS_f2 <- (beta2_f_p2*CS_f2*(alpha2_f_p2))    -(gamma4_p2*ISS_f2)-(omega4_r_f_p2*nu4_f_p2*ISS_f2)+(epsilon2_p2*ISS_f2)-(omega4_d_f_p2*zeta4_f_p2*ISS_f2)+(INF_ISS_f_p2)
  dRR_f2 <- (omega1_r_f_p2*nu1_f_p2*IMR_f2)+(omega2_r_f_p2*nu2_f_p2*ISR_f2)
  dRS_f2 <- (omega3_r_f_p2*nu3_f_p2*IMS_f2)+(omega4_r_f_p2*nu4_f_p2*ISS_f2)
  dDR_f2 <- (omega1_d_f_p2*zeta1_f_p2*IMR_f2)+(omega2_d_f_p2*zeta2_f_p2*ISR_f2)
  dDS_f2 <- (omega3_d_f_p2*zeta3_f_p2*IMS_f2)+(omega4_d_f_p2*zeta4_f_p2*ISS_f2)
  dN_to2<- dU_m2+ dCR_m2+ dCS_m2+ dIMR_m2+ dISR_m2+ dIMS_m2+ dISS_m2 +dU_f2+ dCR_f2+ dCS_f2+ dIMR_f2+ dISR_f2+ dIMS_f2+ dISS_f2
  dutility <- u_healthy*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +u_icu*(ISR_m2+ISS_m2+ISR_f2+ISS_f2)+ u_gw*(IMS_f2+ IMR_f2+IMS_m2+ IMR_m2) + u_healthy*(RR_f2+RR_m2+RS_f2+RS_m2)      
  dcost <-  (Influx_women)*(c_isolation)+c_general_ward*(U_m2+ CR_m2+ CS_m2+U_f2+ CR_f2+ CS_f2) +c_intermediate_ward*(IMR_m2+ IMS_m2+IMR_f2+ IMS_f2)+ c_icu_ward*(ISR_m2+ ISS_m2+ISR_f2+ ISS_f2) 
  dnew_admin <- influx_nonARB + influx_ARB
  ddischarge <- U_m2*Disch_U_m_p2+CR_m2*Disch_CR_m_p2+CS_m2*Disch_CS_m_p2+U_f2*Disch_U_f_p2+CR_f2*Disch_CR_f_p2+CS_f2*Disch_CS_f_p2
  #discharge<- state["discharge"] #list results ddischarge
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2, dutility, dcost, dnew_admin, ddischarge))
  return(results2)
}

# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- SOLVING EQUATIONS--- ------ --- --- ------ --- --- --- --- --- #
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#

#TEST_odes#######
#ode(y = state2, times = times, func = ARB_model_2ch_do_nothing, parms = parameters2, method = "rk4")
#ode(y = state2, times = times, func = ARB_model_2ch_td_newadm, parms = parameters2, method = "rk4")
#ode(y = state2, times = times, func = ARB_model_2ch2_td_newadm, parms = parameters2, method = "rk4")
#ode(y = state2, times = times, func = ARB_model_2pcr_td_newadm, parms = parameters2, method = "rk4")

#ode(y = state2, times = times, func = ARB_model_2ch_tiso_newadm, parms = parameters2, method = "rk4")
#ode(y = state2, times = times, func = ARB_model_2ch2_tiso_newadm, parms = parameters2, method = "rk4")
#ode(y = state2, times = times, func = ARB_model_2pcr_tiso_newadm, parms = parameters2, method = "rk4")

#ode(y = state2, times = times, func = ARB_model_2ch_td_newadmHR_m, parms = parameters2, method = "rk4")
#ode(y = state2, times = times, func = ARB_model_2ch2_td_newadmHR_m, parms = parameters2, method = "rk4")
#ode(y = state2, times = times, func = ARB_model_2_pcr_td_newadmHR_m, parms = parameters2, method = "rk4")

#ode(y = state2, times = times, func = ARB_model_2ch_td_newadmHR_f, parms = parameters2, method = "rk4")
#ode(y = state2, times = times, func = ARB_model_2ch2_td_newadmHR_f, parms = parameters2, method = "rk4")
#ode(y = state2, times = times, func = ARB_model_2_pcr_td_newadmHR_f, parms = parameters2, method = "rk4")

#ode(y = state2, times = times, func = ARB_model_2preE_newadm_all , parms = parameters2, method = "rk4")
#ode(y = state2, times = times, func = ARB_model_2preE_newadm_m, parms = parameters2, method = "rk4")
#ode(y = state2, times = times, func = ARB_model_2preE_newadm_f, parms = parameters2, method = "rk4")


#####

# Define the models you want to run according to the functions (per strategy)
model_list <- c("ARB_model_2ch_do_nothing", "ARB_model_2ch_td_newadm","ARB_model_2ch2_td_newadm","ARB_model_2pcr_td_newadm","ARB_model_2ch_tiso_newadm","ARB_model_2ch2_tiso_newadm","ARB_model_2pcr_tiso_newadm",
                "ARB_model_2ch_td_newadmHR_m","ARB_model_2ch2_td_newadmHR_m","ARB_model_2_pcr_td_newadmHR_m","ARB_model_2ch_td_newadmHR_f","ARB_model_2ch2_td_newadmHR_f","ARB_model_2_pcr_td_newadmHR_f","ARB_model_2preE_newadm_all","ARB_model_2preE_newadm_m","ARB_model_2preE_newadm_f")  # Add last interventions with high-risk groups and pre-emptive isolation
results_epi_ac <- matrix(nrow = 10, ncol = length(model_list))  # Adjust the number of columns based on the number of models/ rows for number of parameters extracted
results_epi_CREprev <- matrix(nrow = 366, ncol = length(model_list))
results_epi_CREinfe <- matrix(nrow = 366, ncol = length(model_list))
results_epi_CREdead <- matrix(nrow = 366, ncol = length(model_list))
results_epi_CREinfe_all <- matrix(nrow = 366, ncol = length(model_list))
results_epi_CREdead_all <- matrix(nrow = 366, ncol = length(model_list))
results_epi_uncolonised <- matrix(nrow = 366, ncol = length(model_list))

results_econ <- matrix(nrow = 3, ncol = length(model_list))  # Adjust the number of columns based on the number of models/ rows for number of parameters extracted
for (i in 1:length(model_list)) {
  #Run initial conditions first:
  times <- seq(from=0, to=365, by = 1)  # Simulate over a year
  #Enterobacterales[CRE/CSE] states####### 
  # ----------------------------------#
  N <- 1000  # Total population size
  # Initial conditions (population sizes in each group)
  U_m20 <- 0.44 * N *(1-0.52)
  CR_m20 <- 0.1445 * N *(1-0.52)
  CS_m20 <-  0.4155* N *(1-0.52)
  IMR_m20 <- 0.09 * CR_m20* (1-0.4283)
  ISR_m20 <- 0.09 * CR_m20* 0.4283
  IMS_m20 <- 0.04 * CS_m20*(1-0.3548)
  ISS_m20 <- 0.04 * CS_m20*0.3548
  RR_m20 <-0
  RS_m20  <-0
  DR_m20  <-0
  DS_m20 <-0
  N_to2<-1050
  utility_to0<-0
  cost_to0<-0
  new_admin0<-0
  discharge0<-0
  
  U_f20 <- 0.44 * N *0.52
  CR_f20 <- 0.1445 * N *0.52
  CS_f20 <-  0.4155* N *0.52
  IMR_f20 <- 0.09 * CR_f20* (1-0.4538)
  ISR_f20 <- 0.09 * CR_f20*0.4538
  IMS_f20 <- 0.04 * CS_f20*(1-0.3832)
  ISS_f20 <- 0.04 * CS_f20*0.3832
  RR_f20 <-0
  RS_f20  <-0
  DR_f20  <-0
  DS_f20 <-0
  
  N_0m20<-  U_m20 + CR_m20 + CS_m20 + IMR_m20 + ISR_m20 + IMS_m20 + ISS_m20 +  RR_m20 + RS_m20 + DR_m20 +DS_m20
  N_0f20<-  U_f20 + CR_f20 + CS_f20 + IMR_f20 + ISR_f20 + IMS_f20 + ISS_f20 +  RR_f20 + RS_f20 + DR_f20 +DS_f20
  N_to0<- N_0m20 + N_0f20
  
  state2 <- c(U_m2 = U_m20, CR_m2=CR_m20, CS_m2= CS_m20, IMR_m2= IMR_m20, ISR_m2=ISR_m20, IMS_m2= IMS_m20, ISS_m2= ISS_m20, RR_m2= RR_m20, RS_m2=RS_m20, DR_m2= DR_m20, DS_m2=DS_m20,
              U_f2 = U_f20, CR_f2=CR_f20, CS_f2= CS_f20, IMR_f2= IMR_f20, ISR_f2=ISR_f20, IMS_f2= IMS_f20, ISS_f2= ISS_f20, RR_f2= RR_f20, RS_f2=RS_f20, DR_f2= DR_f20, DS_f2=DS_f20, N_to2=N_to0, utility=utility_to0, cost=cost_to0, new_admin=new_admin0, discharge=discharge0)
  N_orig2<-N_0m20 + N_0f20
  N_tdif <- N_orig2
  Nt2_spec2<-1030
  #Enterobacterales[CRE/CSE]params####### 
  # ----------------------------------------#
  
  zeta3_m_p2 = 0.228  # mortality rate from IMS  , male
  zeta3_f_p2 = zeta3_m_p2*0.81  # mortality rate from IMS  , female
  zeta1_m_p2 = zeta3_m_p2*1.80  # mortality rate from IMR, male
  zeta1_f_p2 = zeta3_m_p2*0.55  # mortality rate from IMR, female
  zeta2_m_p2 = zeta3_m_p2*1.30  # mortality rate from ISR, male
  zeta2_f_p2 = zeta3_m_p2*2.40  # mortality rate from ISR, female
  zeta4_m_p2 = zeta3_m_p2*1.62  # mortality rate from ISS , male
  zeta4_f_p2 = zeta3_m_p2*2.23  # mortality rate from ISS , female
  kappa_p2 = (35.9-6.2)/(49.5-4.85)
  parameters2 <- c(
    #clearance -natural- parameters
    delta1_p2 = 0.001,  # Clearance value from CR state
    delta2_p2 = 0.001,  # Clearance value from CS state
    #Dicharge rates from uncolonised and colonised
    Disch_U_f_p2 = 1/6,
    Disch_U_m_p2 = 1/6,
    Disch_CR_f_p2 = 1/6,
    Disch_CR_m_p2 = 1/6,
    Disch_CS_f_p2 = 1/6,
    Disch_CS_m_p2 = 1/6,
    #Proportion of women among specific populations (i.e., U, CR, CS, IMR, ISR, IMS, and ISS) [1/unit time] [%]
    mu0_p2 = 0.52,  # % of women among U
    mu1_p2 = 0.52,  # % of women among CR
    mu2_p2 = 0.52,  # % of women among CS
    mu3_p2 = 0.3394,  # % of women among IMR
    mu4_p2 = 0.3280,  # % of women among ISR
    mu5_p2 = 0.50,  # % of women among IMS
    mu6_p2 = 0.44,  # % of women among IMS
    
    #Exposure to anbiotics 
    psi_m_p2 = 0.2225, # % of  individuals exposed  to vancomycin/penicillin among males
    psi_w_p2 = 0.2026,  # % of  individuals exposed to vancomycin/penicillin among males
    
    #Percentage of people under treatment for CRE decolonisation
    psi_mtr_p2=0.2225/1.02,
    psi_wtr_p2=0.2026/1.02,
    
    #Fitness cost. c reduces the transmission rate among resistant strains [1/unit time] [%].
    c_p2=(1-0.927),
    
    #Progression to the development of infection from colonisation among CR and CS states. 
    beta1_m_p2 = (1/22)*0.213, # inverse of LOS plus progression from colonisation to infection among males CR
    beta2_m_p2 = (1/20)*0.034,  # inverse of LOS plus progression from colonisation to infection among males CS
    beta1_f_p2 = (1/27)*0.213, # inverse of LOS plus progression from colonisation to infection among females CR
    beta2_f_p2 = (1/17)*0.034,  # inverse of LOS plus progression from colonisation to infection among females CS
    
    #Natural clearance of mild and severe infections among CR and CS states, respectively [1/unit time] [%].
    gamma1_p2 = 0.001,  # Natural clearance among mild infections R
    gamma2_p2 = 0.001,  # Natural clearance among severe infections R
    gamma3_p2 = 0.001,  # Natural clearance among mild infections S
    gamma4_p2 = 0.001, # Natural clearance among severe infections S
    
    #Mean time of infection considering length of hospital stays [1/length of hospital stay]. 
    omega1_d_m_p2 =(1/21), #IMR patients who died, male
    omega1_r_m_p2 =(1/26), #IMR patients who recovered, male
    omega1_d_f_p2 =(1/31), #IMR patients who died, female
    omega1_r_f_p2 =(1/30), #IMR patients who recovered, female
    omega2_d_m_p2 =(1/7), #ISR patients who died, male
    omega2_r_m_p2 =(1/20),  #ISR patients who recovered, male
    omega2_d_f_p2 = (1/20), #ISR patients who died, female 
    omega2_r_f_p2 =(1/23), #ISR patients who recovered, female
    omega3_d_m_p2 =(1/12), #IMS patients who died, male
    omega3_r_m_p2 =(1/20), #IMS patients who recovered, male
    omega3_d_f_p2 =(1/10), #IMS patients who died, female
    omega3_r_f_p2 =(1/18), #IMS patients who recovered, female
    omega4_d_m_p2 =(1/11), #ISS patients who died, male
    omega4_r_m_p2 =(1/14), #ISS patients who recovered, male
    omega4_d_f_p2 =(1/9), #ISS patients who died, female
    omega4_r_f_p2 =(1/15), #ISS patients who recovered, female
    
    #Percentage of inpatients with CR or CS, respectively, progressing to severe infection in intensive care units [1/unit time] [%].
    alpha1_m_p2 = 0.4283, #% patients with CR progressing to severe infection, males
    alpha2_m_p2 = 0.3548, #% patients with CS progressing to severe infection, males
    alpha1_f_p2 = 0.4585, #% patients with CR progressing to severe infection, females
    alpha2_f_p2 = 0.3832, #% patients with CS progressing to severe infection, females
    
    #Progression from mild to severe infection from IMR and IMS, respectively [1/unit time] [%].
    epsilon1_p2 = 0.01, #progression from IMR to ISR
    epsilon2_p2 = 0.01, #progression from IMS to ISS
    
    #Mortality rates from infection. ζ1 and ζ2 are mortality rates from mild and severe resistant infections, respectively. ζ3 and ζ4 are from mild and severe susceptible infections, respectively [1/unit time] [%].
    zeta3_m_p2 = 0.228,  # mortality rate from IMS  , male
    zeta3_f_p2 = zeta3_m_p2*0.81,  # mortality rate from IMS  , female
    zeta1_m_p2 = zeta3_m_p2*1.80,  # mortality rate from IMR, male
    zeta1_f_p2 = zeta3_m_p2*0.55,  # mortality rate from IMR, female
    zeta2_m_p2 = zeta3_m_p2*1.30,  # mortality rate from ISR, male
    zeta2_f_p2 = zeta3_m_p2*2.40,  # mortality rate from ISR, female
    zeta4_m_p2 = zeta3_m_p2*1.62,  # mortality rate from ISS , male
    zeta4_f_p2 = zeta3_m_p2*2.23,  # mortality rate from ISS , female
    
    #Recovery rates from infection, including IMR, ISR, IMS and ISS due to treatment received [1/unit time] [%].
    nu1_m_p2 = (1-zeta1_m_p2), #recovery rates from IMR, males
    nu1_f_p2 = (1-zeta1_f_p2), #recovery rates from IMR, females
    nu2_m_p2 = (1-zeta2_m_p2), #recovery rates from ISR, males
    nu2_f_p2 = (1-zeta2_f_p2), #recovery rates from ISR, females
    nu3_m_p2 = (1-zeta3_m_p2), #recovery rates from IMS, males
    nu3_f_p2 = (1-zeta3_f_p2), #recovery rates from IMS, females
    nu4_m_p2 =  (1-zeta4_m_p2), #recovery rates from ISS, males
    nu4_f_p2 =  (1-zeta4_f_p2),#recovery rates from ISS, females
    
    #Constant background rate that captures transmission from non-human sources, horizontal transmission, or de novo emergence [1/unit time] [number].
    b_p2 = 0.01,
    
    # Percentage of people with resistant infections receiving inappropriate empirical antibiotic treatment [1/unit time] [%].
    phi_m_p2 = 0.3782,  # Placeholder value, adjust as needed
    phi_f_p2 = 0.3846,  # Placeholder value, adjust as needed
    
    #Factor of burden associated with inappropriate empirical antibiotic treatment and increased ICU admission among resistant infections [1/unit time] [%].
    pi_p2= 1.02,
    
    #Transmission parameter {update this correspondingly after calibrating it with real data}
    tau_p2= 0.3986551,
    
    #community-acquired infection upon hospital admission rate
    caIha_p2=0.007,
    
    #percentage of people tested
    test_p2=0.20, 
    HR_perc2=0.2,
    or_HR_scenar1_a=1.04,
    or_HR_scenarMen_a=2.27,
    ### ### ### ### ###
    #sensitivity chrom_1
    sens_chrom_a=0.826,
    #sensitivity chrom_1
    sens_chrom2_a=0.90, 
    #sensitivity chrom_1
    sens_pcr_a=1,
    #turnaround chrom_1
    turn_chrom_a=3,
    #turnaround chrom_1
    turn_chrom2_a=2,
    #turnaround pcr_1
    turn_pcr_a=1,
    #isolation contact precaution transmission reduction
    reduc_conpre_a=0.35,
    #efficiency decolonisation
    eff_decol_a= 0.26,
    #effect on self-infection decolonisation
    eff_decol_selfi_a=0.041,
    #Turnaround decolonisation program in days
    turnaround_decol_a=7,
    ##
    #costs wards
    c_general_ward= 50,
    c_intermediate_ward=92,
    c_icu_ward=218,
    c_decol_1pd=72.88,
    c_isolation=42.3,
    c_chrom=10.2,
    c_chrom2=13.6, 
    c_pcr=33,
    c_bc=16.9,
    #utilities
    u_healthy=0.92,
    u_icu=0.92-0.34,
    u_gw=0.64,
    u_recovICU=0.74
    
  )
  
  #####
  # Construct the function name from the model_list string
  model_func_name <- get(model_list[i])
  # Solve the ODE using the 'deSolve' package's ode function
  O_solution <- ode(y = state2, times = times, func = model_func_name, parms = parameters2, method = "rk4")
  O_solution2 <- as.data.frame(O_solution)
  O_solution <- as.data.frame(O_solution)
  # Process the ODE solution: sum all values reported in the first 350 rows for your variable of interest
  results_epi_ac[1, i] <- sum(O_solution[1:366, "CR_m2"]) + sum(O_solution[1:366, "CR_f2"]) + sum(O_solution[1:366, "IMR_m2"]) + sum(O_solution[1:366, "IMR_f2"])+ sum(O_solution[1:366, "ISR_m2"])+sum(O_solution[1:366, "ISR_f2"])
  results_epi_ac[2, i] <- sum(O_solution[1:366, "IMR_m2"]) + sum(O_solution[1:366, "IMR_f2"])+ sum(O_solution[1:366, "ISR_m2"])+sum(O_solution[1:366, "ISR_f2"])
  results_epi_ac[3, i] <- (O_solution[366, "DR_f2"]) + (O_solution[366, "DR_m2"]) 
  results_epi_ac[4, i] <- sum(O_solution[1:366, "IMR_m2"]) + sum(O_solution[1:366, "IMR_f2"])+ sum(O_solution[1:366, "ISR_m2"])+sum(O_solution[1:366, "ISR_f2"]) + sum(O_solution[1:366, "IMS_m2"]) + sum(O_solution[1:366, "IMS_f2"])+ sum(O_solution[1:366, "ISS_m2"])+sum(O_solution[1:366, "ISS_f2"])
  results_epi_ac[5, i] <- (O_solution[366, "DR_f2"]) + (O_solution[366, "DR_m2"]) + (O_solution[366, "DS_f2"]) + (O_solution[366, "DS_m2"]) 
  results_epi_ac[6, i] <- sum(O_solution[366, "new_admin"]) 
  results_epi_ac[7, i] <- sum(O_solution[1:366, "U_m2"]) +sum(O_solution[1:366, "U_f2"]) 
  results_epi_ac[8, i] <- sum(O_solution[1:366, "U_m2"]) +sum(O_solution[1:366, "U_f2"]) +sum(O_solution[1:366, "CR_m2"]) +sum(O_solution[1:366, "CR_f2"]) +sum(O_solution[1:366, "CS_m2"]) +sum(O_solution[1:366, "CS_f2"]) 
  results_epi_ac[9, i] <- O_solution[366, "RR_m2"] + O_solution[366, "RR_f2"]+O_solution[366, "RS_m2"] + O_solution[366, "RS_f2"]
  results_epi_ac[10, i] <- (O_solution[366, "discharge"]) 
  
  results_epi_CREprev[, i] <- O_solution2$CR_m2 + O_solution2$CR_f2 + O_solution2$IMR_m2 + O_solution2$IMR_f2 + O_solution2$ISR_m2 + O_solution2$ISR_f2
  results_epi_CREinfe[, i]  <- O_solution2$IMR_m2 + O_solution2$IMR_f2 + O_solution2$ISR_m2 + O_solution2$ISR_f2
  results_epi_CREdead[, i]  <- O_solution2$DR_m2 + O_solution2$DR_f2
  results_epi_CREinfe_all[, i]  <-  O_solution2$IMR_m2 + O_solution2$IMR_f2 + O_solution2$ISR_m2 + O_solution2$ISR_f2 + O_solution2$IMS_m2 + O_solution2$IMS_f2 + O_solution2$ISS_m2 + O_solution2$ISS_f2
  results_epi_CREdead_all[, i]  <- O_solution2$DR_m2 + O_solution2$DR_f2 + O_solution2$DS_m2 + O_solution2$DS_f2
  
  #Economics
  #Store econ results per strategy
  results_econ[1, i] <- O_solution[366, "cost"]
  results_econ[2, i] <- (results_epi_ac[8, i]*0.92)+(sum(O_solution[1:366, "ISR_m2"]))*0.58+(sum(O_solution[1:366, "ISR_f2"])*0.58)+((sum(O_solution[1:366, "IMR_m2"])+sum(O_solution[1:366, "IMR_f2"]))*0.64)+((sum(O_solution[1:366, "IMS_m2"])+sum(O_solution[1:366, "IMS_f2"]))*0.64)+(sum(O_solution[1:366, "ISS_m2"]))*0.58+(sum(O_solution[1:366, "ISS_f2"])*0.58)+ (O_solution[366, "RR_m2"])*0.92+(O_solution[366, "RR_f2"]*0.92)+ (O_solution[366, "RS_m2"])*0.92+(O_solution[366, "RS_f2"]*0.92)+results_epi_ac[10, i]*0.92
  results_econ[3, i] <- results_econ[2, i]/O_solution[366,"new_admin"] #check if usage is appropriate
}
#Compute ICER per strategy
results_icer <- matrix(nrow = 13, ncol = length(model_list))
results_icer[1,1] <- 0
results_icer[2,1] <- results_econ[1, 1]
results_icer[3,1] <- results_econ[2, 1]
results_icer[4,1] <- results_epi_ac[1, 1]
results_icer[5,1] <- results_epi_ac[2, 1]
results_icer[6,1] <- results_epi_ac[3, 1]
results_icer[7,1] <- results_epi_ac[4, 1]
results_icer[8,1] <- results_epi_ac[5, 1]
results_icer[9,1] <- results_epi_ac[9, 1]
results_icer[10,1] <- results_epi_ac[8, 1]
results_icer[11,1] <- results_epi_ac[10, 1]
results_icer[12,1] <- results_epi_ac[10, 1]+  results_epi_ac[9, 1]+results_epi_ac[8, 1]+results_epi_ac[4, 1]
results_icer[13,1] <- results_epi_ac[6, 1]
for (i in 2:length(model_list)) {
  results_icer[1, i] <- 0
  results_icer[2, i] <-  results_econ[1, i]
  results_icer[3, i] <-  results_econ[2, i]
  results_icer[4, i]<- results_epi_ac[1, i]
  results_icer[5, i]<- results_epi_ac[2, i]
  results_icer[6, i]<- results_epi_ac[3, i]
  results_icer[7, i]<- results_epi_ac[4, i]
  results_icer[8, i]<- results_epi_ac[5, i]
  results_icer[9, i]<- results_epi_ac[9, i]
  results_icer[10, i]<- results_epi_ac[8, i]
  results_icer[11, i]<- results_epi_ac[10, i]
  results_icer[12, i]<- results_epi_ac[10, i]+  results_epi_ac[9, i]+results_epi_ac[8, i]+results_epi_ac[4, i]
  results_icer[3, i] <- results_econ[2, i] + ifelse((results_icer[12, 1] - results_icer[12, i]) > 0, (results_icer[12, 1] - results_icer[12, i]) * 0.92, 0)
  results_icer[1, i] <- (results_icer[2, i]-results_icer[2, 1])/(results_icer[3, i]-results_icer[3, 1])
  results_icer[13,i] <- results_epi_ac[6, i]
  
}
results_epi_acC<- results_epi_ac
results_econC<- results_econ
results_icerC<- results_icer
transposed_icer <- t(results_icer)
transposed_icerdf <- as.data.frame(transposed_icer)
colnames(transposed_icerdf) <- c("ICER", "Costs", "QALYs","ARB colonisation","ARB infections", "ARB deaths", "Total infections","Total deaths","Total Recovered","Total U,CR,CS", "Total discharge","total population")  # Add more names as needed
rownames(transposed_icerdf) <- c("Do-nothing", "Strategy2", "Strategy3", "Strategy4", "Strategy5", "Strategy6", "Strategy7", "Strategy8", "Strategy9", "Strategy10","Strategy11", "Strategy12","Strategy13","Strategy14","Strategy15", "Strategy16")
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/AMR_BSI_BurdenChile/0_Article_AMR Transmission dynamics Chile & LMICs/0_BSIModelling/0_analysis/0_Figures_model")
write.csv(transposed_icerdf, "matrix_resCRE.csv", row.names = TRUE)

colonCRE_a <- matrix(nrow =1, ncol =  length(model_list))
infecCRE_a <- matrix(nrow =1, ncol =  length(model_list))
DeadCRE_a <- matrix(nrow =1, ncol =  length(model_list))
Dead_CRE_tot <- matrix(nrow =1, ncol =  length(model_list))
for (i in 1:length(model_list)) {
  colonCRE_a[1, i] <- ((results_epi_acC[1, i]-results_epi_acC[2, i])/results_epi_acC[6, i])*100
  infecCRE_a[1, i] <- (results_epi_acC[2, i]/results_epi_acC[6, i])*100
  DeadCRE_a[1, i] <-  (results_epi_acC[3, i]/results_epi_acC[6, i])*100
  Dead_CRE_tot[1, i] <-  (results_epi_ac[3, i])
}
# Convert matrices to vectors
colonCRE_a <- as.vector(t(colonCRE_a))  # Transpose and then convert to vector
infecCRE_a <- as.vector(t(infecCRE_a))  # Assuming similar structure
DeadCRE_a<- as.vector(t(DeadCRE_a))    # Assuming similar structure
Dead_CRE_tot<- as.vector(t(Dead_MRSA_tot))

