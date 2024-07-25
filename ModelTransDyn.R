# Load required library
library(ggplot2)
library(deSolve)
library(ggplot2)
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/AMR_BSI_BurdenChile/0_Article_AMR Transmission dynamics Chile & LMICs/0_BSIModelling/0_analysis")

##################################################################################################
##################################################################################################
######MRSA#######
##################################################################################################
##################################################################################################

# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- MODEL & EQUATIONS BELOW --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# Define a function for the differential equations
ARB_model_1 <- function(times, state, parms){
  
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
  test_p1 <-parms["test"]
  
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
  
  results1 <- list(c(dU_m1, dCR_m1, dCS_m1, dIMR_m1, dISR_m1, dIMS_m1, dISS_m1, dRR_m1, dRS_m1, dDR_m1, dDS_m1,
                     dU_f1, dCR_f1, dCS_f1, dIMR_f1, dISR_f1, dIMS_f1, dISS_f1, dRR_f1, dRS_f1, dDR_f1, dDS_f1, dN_to))
  
  return(results1)
  
}
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
  #percentage of people tested
  test_p1=0
)




# ----------------------------------------#
######
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- ---   TIME  --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# Define the time span for simulation
times <- seq(from=0, to=365, by = 1)  # Simulate over a year

# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- ---BASELINE CONDITIONS BELOW--- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# ----------------------------------------#
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
            U_f1 = U_f10, CR_f1=CR_f10, CS_f1= CS_f10, IMR_f1= IMR_f10, ISR_f1=ISR_f10, IMS_f1= IMS_f10, ISS_f1= ISS_f10, RR_f1= RR_f10, RS_f1=RS_f10, DR_f1= DR_f10, DS_f1=DS_f10, N_to=N_to0)

N_orig1<-N_0m10 + N_0f10
N_tdif <- N_orig1

# ----------------------------------#
#####

# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- SOLVING EQUATIONS--- ------ --- --- ------ --- --- --- --- --- #
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# Solving differential equations
library(reshape)
library(ggplot2)

# --------------- # 
#      MRSA
# --------------- #
mrsa_model_output <- ode(y = state1, times = times, func = ARB_model_1, parms =  parameters1,
                         method = "rk4")
mrsa_df <- as.data.frame(mrsa_model_output)
mrsa_df$N_t <- mrsa_df$U_m1+mrsa_df$U_f1+mrsa_df$CR_m1+mrsa_df$CR_f1+mrsa_df$CS_m1+mrsa_df$CS_f1
+mrsa_df$IMR_m1+mrsa_df$IMR_f1+mrsa_df$IMS_m1+mrsa_df$IMS_f1+mrsa_df$ISR_m1+mrsa_df$ISR_f1+
  +mrsa_df$ISS_m1+mrsa_df$ISS_f1
#mrsa_df_long_f <- mrsa_df_long[!mrsa_df_long$variable %in% c("RS_m1", "RR_m1", "RS_f1", "RR_f1"), ]
mrsa_df$U <- mrsa_df$U_m1 + mrsa_df$U_f1
mrsa_df$CR <- mrsa_df$CR_m1 + mrsa_df$CR_f1
mrsa_df$CS <- mrsa_df$CS_m1 + mrsa_df$CS_f1
mrsa_df$IMR <- mrsa_df$IMR_m1 + mrsa_df$IMR_f1
mrsa_df$ISR <- mrsa_df$ISR_m1 + mrsa_df$ISR_f1
mrsa_df$IMS <- mrsa_df$IMS_m1 + mrsa_df$IMS_f1
mrsa_df$ISS <- mrsa_df$ISS_m1 + mrsa_df$ISS_f1
mrsa_df$RR <- mrsa_df$RR_m1 + mrsa_df$RR_f1
mrsa_df$RS <- mrsa_df$RS_m1 + mrsa_df$RS_f1
mrsa_df$DR <- mrsa_df$DR_m1 + mrsa_df$DR_f1
mrsa_df$DS <- mrsa_df$DS_m1 + mrsa_df$DS_f1
# Filtering for time <= 25 and excluding RS and RR totals for the plot
mrsa_df_f <- mrsa_df[mrsa_df$time <= 350, ]
mrsa_df_long_f <- melt(mrsa_df_f, id.vars = "time", measure.vars = c("U","CR", "CS", "IMR", "ISR", "IMS", "ISS", "DR", "DS"))

# Plotting pop dynamics
tiff("EquilibriumDynamics_MRSA.tiff", width = 13, height = 10, units = 'in', res = 500)
ggplot(mrsa_df_long_f, aes(x = time, y = value, colour = variable)) + 
  geom_line(linewidth = 1.2)+
  theme_minimal() +
  labs(x = "Time (days)", 
       y = "Population size", 
       title = "Dynamics of sex-based MRSA model compartments over time",
       colour = "Compartment") +
  scale_x_continuous(breaks = seq(0, 350, by = 20))+
  scale_y_continuous(breaks = seq(0, 1000, by = 100), limits = c(0, 1000)) +
  theme(legend.position = c(0.9, 1), # Adjust these values to move the legend inside the plot
        legend.justification = c(0, 1))+
  theme(
    axis.text.y = element_text(size = 12, family = "Times New Roman"),  # Set font family for y-axis text
    axis.text.x = element_text(size = 12,family = "Times New Roman"),  # Set font family for x-axis text
    axis.title.y = element_text(family = "Times New Roman", margin = margin(r = 20)),  # Set font family for y-axis title
    axis.title.x = element_text(family = "Times New Roman"),  # Set font family for x-axis title
    text = element_text(size=15, family = "Times New Roman")  # Set font family for all text
  )# Anchors the legend position; (0,1) is top-left
dev.off()


tiff("Popsize_MRSA.tiff", width = 13, height = 10, units = 'in', res = 500)
#Plotting population size over time
ggplot(data = mrsa_df, aes(x = time, y = N_t)) + 
  geom_line(linewidth = 1.5)+
  geom_point() + # Add points
  theme_minimal() + # Use a minimal theme
  labs(title = "Total hospital population size over time",
       x = "Time (days)",
       y = "Total population size (Nt), MRSA") +
  theme(plot.title = element_text(hjust = 0.5))+ # Center the plot title
  scale_y_continuous(breaks = seq(0, 1100, by = 50), limits = c(800, 1100))+
  theme(
    axis.text.y = element_text(size = 12, family = "Times New Roman"),  # Set font family for y-axis text
    axis.text.x = element_text(size = 12,family = "Times New Roman"),  # Set font family for x-axis text
    axis.title.y = element_text(family = "Times New Roman", margin = margin(r = 20)),  # Set font family for y-axis title
    axis.title.x = element_text(family = "Times New Roman"),  # Set font family for x-axis title
    text = element_text(size = 12, family = "Times New Roman")  # Set font family for all text
  )
dev.off()



##################################################################################################
##################################################################################################
######CRE#######
##################################################################################################
##################################################################################################

# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- MODEL & EQUATIONS BELOW --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# Define a function for the differential equations
#N original baseline conditions
ARB_model_2 <- function(times, state, parms) {
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
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_m2)+(psi_wtr_p2*CR_f2)
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
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2))
  return(results2)
}
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- ---  MODEL PARAMETERS BELOW --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# ----------------------------------------#
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
  
  test_p2=0
  
)










######
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- ---   TIME  --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# Define the time span for simulation
times <- seq(from=0, to=365, by = 1)  # Simulate over a year
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- ---BASELINE CONDITIONS BELOW--- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
#####
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
            U_f2 = U_f20, CR_f2=CR_f20, CS_f2= CS_f20, IMR_f2= IMR_f20, ISR_f2=ISR_f20, IMS_f2= IMS_f20, ISS_f2= ISS_f20, RR_f2= RR_f20, RS_f2=RS_f20, DR_f2= DR_f20, DS_f2=DS_f20, N_to2=N_to0)
N_orig2<-N_0m20 + N_0f20
N_tdif <- N_orig2
Nt2_spec2<-1030
#####

# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- SOLVING EQUATIONS--- ------ --- --- ------ --- --- --- --- --- #
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --------------- # 
#      CRE
# --------------- #
# Solving differential equations
library(reshape)
library(ggplot2)
cre_model_output <- ode(y = state2, times = times, func = ARB_model_2, parms =  parameters2,
                        method = "rk4")
cre_df <- as.data.frame(cre_model_output)
cre_df$N_t <- cre_df$U_m2+cre_df$U_f2+cre_df$CR_m2+cre_df$CR_f2+cre_df$CS_m2+cre_df$CS_f2+cre_df$IMR_m2+cre_df$IMR_f2+cre_df$IMS_m2+cre_df$IMS_f2+cre_df$ISR_m2+cre_df$ISR_f2+cre_df$ISS_m2+cre_df$ISS_f2
#mrsa_df_long_f <- mrsa_df_long[!mrsa_df_long$variable %in% c("RS_m1", "RR_m1", "RS_f1", "RR_f1"), ]
cre_df$U <- cre_df$U_m2 + cre_df$U_f2
cre_df$CR <- cre_df$CR_m2 + cre_df$CR_f2
cre_df$CS <- cre_df$CS_m2 + cre_df$CS_f2
cre_df$IMR <- cre_df$IMR_m2 + cre_df$IMR_f2
cre_df$ISR <- cre_df$ISR_m2 + cre_df$ISR_f2
cre_df$IMS <- cre_df$IMS_m2 + cre_df$IMS_f2
cre_df$ISS <- cre_df$ISS_m2 + cre_df$ISS_f2
cre_df$RR <- cre_df$RR_m2 + cre_df$RR_f2
cre_df$RS <- cre_df$RS_m2 + cre_df$RS_f2
cre_df$DR <- cre_df$DR_m2 + cre_df$DR_f2
cre_df$DS <- cre_df$DS_m2 + cre_df$DS_f2
# Filtering for time <= 25 and excluding RS and RR totals for the plot
cre_df_f <- cre_df[cre_df$time <= 350, ]
cre_df_long_f <- melt(cre_df_f, id.vars = "time", measure.vars = c("U","CR", "CS", "IMR", "ISR", "IMS", "ISS", "DR", "DS"))

# Plotting pop dynamics
tiff("EquilibriumDynamics_CRE.tiff", width = 13, height = 10, units = 'in', res = 500)
ggplot(cre_df_long_f, aes(x = time, y = value, colour = variable)) + 
  geom_line(linewidth = 1.2)+
  theme_minimal() +
  labs(x = "Time (days)", 
       y = "Population size", 
       title = "Dynamics of sex-based CRE model compartments over time",
       colour = "Compartment") +
  scale_x_continuous(breaks = seq(0, 350, by = 20))+
  scale_y_continuous(breaks = seq(0, 1000, by = 100), limits = c(0, 1000)) +
  theme(legend.position = c(0.9, 1), # Adjust these values to move the legend inside the plot
        legend.justification = c(0, 1))+
  theme(
    axis.text.y = element_text(size = 12, family = "Times New Roman"),  # Set font family for y-axis text
    axis.text.x = element_text(size = 12,family = "Times New Roman"),  # Set font family for x-axis text
    axis.title.y = element_text(family = "Times New Roman", margin = margin(r = 20)),  # Set font family for y-axis title
    axis.title.x = element_text(family = "Times New Roman"),  # Set font family for x-axis title
    text = element_text(size=15, family = "Times New Roman")  # Set font family for all text
  )# Anchors the legend position; (0,1) is top-left
dev.off()


tiff("Popsize_CRE.tiff", width = 13, height = 10, units = 'in', res = 500)
#Plotting population size over time
ggplot(data = cre_df, aes(x = time, y = N_t)) + 
  geom_line(linewidth = 1.5)+
  geom_point() + # Add points
  theme_minimal() + # Use a minimal theme
  labs(title = "Total hospital population size over time",
       x = "Time (days)",
       y = "Total population size (Nt), CRE") +
  theme(plot.title = element_text(hjust = 0.5))+ # Center the plot title
  scale_y_continuous(breaks = seq(0, 1100, by = 50), limits = c(800, 1100))+
  theme(
    axis.text.y = element_text(size = 12, family = "Times New Roman"),  # Set font family for y-axis text
    axis.text.x = element_text(size = 12,family = "Times New Roman"),  # Set font family for x-axis text
    axis.title.y = element_text(family = "Times New Roman", margin = margin(r = 20)),  # Set font family for y-axis title
    axis.title.x = element_text(family = "Times New Roman"),  # Set font family for x-axis title
    text = element_text(size = 12, family = "Times New Roman")  # Set font family for all text
  )
dev.off()


########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################

# https://cran.r-project.org/web/packages/sensitivity/sensitivity.pdf



########################################################################
########################################################################
########################################################################
#Sensitivity analyses: GLOBAL SENS anALysIS
########################################################################
########################################################################
########################################################################
library(sensitivity)
library(deSolve)
library(lhs)
################################################
#I. MRSA
################################################
# Calculate the min and max for each parameter
param_ranges1 <- sapply(parameters1, function(x) c(x * 0.9, x * 1.1))
# Convert the param_ranges to a matrix for easier handling in LHS
param_ranges_matrix1 <- matrix(unlist(param_ranges1), nrow = length(parameters1), byrow = TRUE, dimnames = list(names(parameters1), c("min", "max")))
params_range1 <- data.frame(t(param_ranges1)) 
# Define the parameter names and their respective ranges for sensitivity analysis
#####
# Number of samples
n_samples <- 100
# Generate LHS samples
lhs_samples1 <- randomLHS(n = n_samples, k = length(parameters1))
# Scale the LHS samples to your parameters' ranges
scaled_samples1 <- apply(lhs_samples1, MARGIN = 2, FUN = function(x, range) {
  range[1] + (range[2] - range[1]) * x
}, range = t(params_range1))

simulation_results1 <- vector("list", n_samples)
times <- seq(from = 0, to = 365, by = 1)  # Adjust time range as needed

for (i in 1:n_samples) {
  parms_i <- setNames(as.list(scaled_samples1[i, ]), names(parameters1))
  parms_i <- unlist(list((parms_i)))
  simulation_results1[[i]] <- ode(y = state1, times = times, func = ARB_model_1, parms = parms_i, method = "rk4")
}

#create outputs from above.
final_values_X1a <- numeric(n_samples)  # Replace 'n_samples' with the actual number of samples
final_values_X2b <- numeric(n_samples)  # Replace 'n_samples' with the actual number of samples
final_values_X3c <- numeric(n_samples)  # Replace 'n_samples' with the actual number of samples

for (i in 1:n_samples) {
  # Assuming 'X' is the name of your state variable and it's stored in a column named 'X'
  # Replace 'X' with the actual name or index of your variable
  final_values_X1a[i] <- simulation_results1[[i]][nrow(simulation_results1[[i]]), "CR_m1"]+simulation_results1[[i]][nrow(simulation_results1[[i]]), "CR_f1"]+ simulation_results1[[i]][nrow(simulation_results1[[i]]), "IMR_m1"]+simulation_results1[[i]][nrow(simulation_results1[[i]]), "IMR_f1"]+ simulation_results1[[i]][nrow(simulation_results1[[i]]), "ISR_m1"]+simulation_results1[[i]][nrow(simulation_results1[[i]]), "ISR_f1"]
  final_values_X2b[i] <- simulation_results1[[i]][nrow(simulation_results1[[i]]), "IMR_m1"]+simulation_results1[[i]][nrow(simulation_results1[[i]]), "IMR_f1"]+ simulation_results1[[i]][nrow(simulation_results1[[i]]), "ISR_m1"]+simulation_results1[[i]][nrow(simulation_results1[[i]]), "ISR_f1"]
  final_values_X3c[i] <- simulation_results1[[i]][nrow(simulation_results1[[i]]), "DR_m1"]+simulation_results1[[i]][nrow(simulation_results1[[i]]), "DR_f1"]
}
results_df_allR1 <- data.frame(X_final_value = final_values_X1a)
results_df_allRinf1 <- data.frame(X_final_value = final_values_X2b)
results_df_deathsR1 <- data.frame(X_final_value = final_values_X3c)

lhs_ranks1 <- apply(scaled_samples1 , 2, rank)
names(lhs_ranks1)<-colnames(param_names1)
output_ranks1 <- apply(results_df_allR1, 2, rank)
output_ranks1 <- as.numeric(output_ranks1)  # Ensuring it's a numeric vector

library(corrr)
# Compute PRCC between each parameter and each output
prcc_results1 <- cor(lhs_ranks1, output_ranks1)
print(prcc_results1)
prcc_values1 <- numeric(ncol(lhs_ranks1))
for (i in 1:ncol(lhs_ranks1)) {
  prcc_values1[i] <- cor(lhs_ranks1[, i], output_ranks1, method = "spearman")
}
# Initialize a vector to hold p-values
prcc_p_values1 <- numeric(ncol(lhs_ranks1))
# Perform the tests
for (i in 1:ncol(lhs_ranks1)) {
  test_result1 <- cor.test(lhs_ranks1[, i], output_ranks1, method = "spearman")
  prcc_p_values1[i] <- test_result1$p.value
}
#ParameterLabelNames####
param_names1 <- c(
  "δ1",  
  "δ2",  
  "Discharge U[w]",
  "Discharge U[m]",
  "Discharge CR[w]",
  "Discharge CR[m]",
  "Discharge CS[w]",
  "Discharge CS[m]",
  "μ0",
  "μ1",
  "μ2",
  "μ3",
  "μ4",
  "μ5",
  "μ6",
  "Ψ[m]",
  "Ψ[w]",
  "Ψt[m]",
  "Ψt[w]",
  "c",
  "β1[m]",
  "β2[m]",
  "β1[w]",
  "β2[w]",
  "γ1",
  "γ2",
  "γ3",
  "γ4",
  "ω1d[m]",
  "ω1r[m]",
  "ω1d[w]",
  "ω1r[w]",
  "ω2d[m]",
  "ω2r[m]",
  "ω2d[w]",
  "ω2r[w]",
  "ω3d[m]",
  "ω3r[m]",
  "ω3d[w]",
  "ω3r[w]",
  "ω4d[m]",
  "ω4r[m]",
  "ω4d[w]",
  "ω4r[w]",
  "α1[m]",
  "α2[m]",
  "α1[w]",
  "α2[w]",
  "ε1",
  "ε2",
  "ζ3[m]",
  "ζ3[w]",
  "ζ1[m]",
  "ζ1[w]",
  "ζ2[m]",
  "ζ2[w]",
  "ζ4[m]",
  "ζ4[w]",
  "ν1[m]",
  "ν1[w]",
  "ν2[m]",
  "ν2[w]",
  "ν3[m]",
  "ν3[w]",
  "ν4[m]",
  "ν4[w]",
  "b",
  "φ[m]",
  "φ[w]",
  "π",
  "τ",
  "caIha", 
  "test"
)
#####
# Create the data frame using the named vectors
prcc_results1 <- data.frame(
  Parameter = param_names1,
  PRCC = prcc_values1,
  P_Value = prcc_p_values1  # Include this line only if you have computed p-values
)
prcc_results1$P_Value <- ifelse(prcc_results1$P_Value < 0.001, "<0.001", sprintf("%.3f", prcc_results1$P_Value))

#Latin Hypercube Sampling/Partial Rank Correlation Coefficient (LHS/PRCC) 
## Sort prcc_results based on the absolute values of PRCC but keep the original sign
prcc_results_sorted1 <- prcc_results1[order(abs(prcc_results1$PRCC), decreasing = TRUE), ]
library(ggplot2)
library(scales)
tiff("GlobalsensA_MRSA.tiff", width = 13, height = 10, units = 'in', res = 500)
ggplot(prcc_results1, aes(x = reorder(Parameter, PRCC), y = PRCC, fill = PRCC)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flips the axes for easier reading
  scale_fill_gradient2(low = "cornflowerblue", mid = "beige", high = "gold2", midpoint = 0) +  # Color gradient
  geom_text(aes(label = paste("p ", P_Value), y = ifelse(PRCC < 0, PRCC - max(abs(PRCC) / 30), PRCC + max(abs(PRCC) / 30))), 
            position = position_dodge(width = 0.9), 
            hjust = ifelse(prcc_results1$PRCC < 0, 1.1, -0.1),  # Adjust text position based on PRCC value
            size = 3, family = "Times New Roman") +  # Adjust text size as needed
  labs(title = "",
       x = "Parameters",
       y = "PRCC value") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(-0.7, 0.6, by = 0.1), labels = scales::label_number(auto = TRUE))+ 
  theme(
    legend.position = "none",  # Remove legend if unnecessary
    axis.text.y = element_text(size = 10.5, family = "Times New Roman"),  # Set font family for y-axis text
    axis.text.x = element_text(size = 10.5, family = "Times New Roman"),  # Set font family for x-axis text
    axis.title.y = element_text(family = "Times New Roman", margin = margin(r = 20)),  # Set font family for y-axis title
    axis.title.x = element_text(family = "Times New Roman"),  # Set font family for x-axis title
    text = element_text(family = "Times New Roman")  # Set font family for all text
  )
dev.off()  # Close the


library(ggplot2)
library(dplyr)
# Filter the results for significant p-values
significant_prcc_results <- prcc_results1 %>%
  filter(P_Value <= 0.05)
# Now plot using the filtered results
library(ggplot2)
# Plot using the filtered results
library(ggplot2)
# Plot using the filtered results
glob_sens_mrsa_pval <- ggplot(significant_prcc_results, aes(x = reorder(Parameter, PRCC), y = PRCC, fill = PRCC)) +
  geom_bar(stat = "identity", color = "black") +  # Add black border around bars
  coord_flip() +  # Flips the axes for easier reading
  scale_fill_gradient2(low = "cornflowerblue", mid = "beige", high = "gold2", midpoint = 0) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) +  # Add horizontal line at y=0 (which is x=0 after flipping)
  geom_text(aes(label = paste("p = ", format(P_Value, digits = 2)), 
                y = ifelse(PRCC < 0, PRCC - max(abs(PRCC) / 30), PRCC + max(abs(PRCC) / 30))), 
            position = position_dodge(width = 0.1), 
            hjust = ifelse(significant_prcc_results$PRCC < 0, 1.1, -0.1), 
            size = 5, family = "Times New Roman") +
  labs(title = "",
       x = "Parameters",
       y = "PRCC value") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 21, family = "Times New Roman"),
        axis.text.x = element_text(size = 21, family = "Times New Roman"),
        axis.title.y = element_text(family = "Times New Roman", size = 24, margin = margin(r = 20)),
        axis.title.x = element_text(family = "Times New Roman", size = 24),
        text = element_text(family = "Times New Roman")) +
  scale_y_continuous(breaks = seq(-0.8, 0.8, by = 0.2), labels = scales::label_number(auto = TRUE), limits = c(-0.9, 0.9))  # Extend y-axis limits
# Save the plot
ggsave("glob_sens_mrsa_pval.tiff", glob_sens_mrsa_pval, width = 12, height = 7, dpi = 800)





################################################
#II. CRE
################################################

ARB_model_2 <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
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
  
  N2_2 <- U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 + RR_f2 + RS_f2 + DR_f2 + DS_f2
  
  #N total (women+men)
  Nt_2 = N1_2 + N2_2 
  #population at time t
  Nt2_spec2 <-  U_f2 + CR_f2 + CS_f2 + IMR_f2 + ISR_f2 + IMS_f2 + ISS_f2 +U_m2 + CR_m2 + CS_m2 + IMR_m2 + ISR_m2 + IMS_m2 + ISS_m2 
  if(!is.numeric(Nt2_spec2)) stop("Nt2_spec2 is not numeric")
  
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
  
  #N original baseline conditions
  N_orig2<-N_0m20 + N_0f20
  if(!is.numeric(Nt2_spec2)) stop("Nt2_spec2 is not numeric")
  if(any(!sapply(parms, is.numeric))) stop("One or more parameters are not numeric")
  
  #Influx of populations
  INF_U_f_p2 <- (1050-Nt2_spec2)*(0.44)*mu0_p2
  INF_U_m_p2 <- (1050-Nt2_spec2)*0.44*(1-mu0_p2)
  INF_CR_f_p2 <- (1050-Nt2_spec2)*0.1445*mu1_p2
  INF_CR_m_p2 <- (1050-Nt2_spec2)*0.1445*(1-mu1_p2)
  INF_CS_f_p2 <- (1050-Nt2_spec2)*0.4155*mu2_p2
  INF_CS_m_p2 <- (1050-Nt2_spec2)*0.4155*(1-mu2_p2)
  INF_IMR_f_p2 <- (1050-Nt2_spec2)*(caIha_p2)*mu3_p2*(1/8)
  INF_IMR_m_p2<- (1050-Nt2_spec2)*(caIha_p2)*(1-mu3_p2)*(1/8)
  INF_ISR_f_p2<- (1050-Nt2_spec2)*(caIha_p2)*mu4_p2*(1/8)
  INF_ISR_m_p2<- (1050-Nt2_spec2)*(caIha_p2)*(1-mu4_p2)*(1/8)
  INF_IMS_f_p2<- (1050-Nt2_spec2)*(caIha_p2)*(mu5_p2)*(1/8)
  INF_IMS_m_p2<- (1050-Nt2_spec2)*(caIha_p2)*(1-mu5_p2)*(1/8)
  INF_ISS_f_p2<- (1050-Nt2_spec2)*(caIha_p2)*(mu6_p2)*(1/8)
  INF_ISS_m_p2<- (1050-Nt2_spec2)*(caIha_p2)*(1-mu6_p2)*(1/8)
  
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
  
  dU_f2 <-   (delta1_p2*CR_f2)+(delta2_p2*CS_f2)+(INF_U_f_p2)-(FOC_u_2*mu0_p2)-(U_f2*Disch_U_f_p2)+(psi_w_p2*CS_m2)+(psi_wtr_p2*CR_f2)
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
  
  results2 <- list(c(dU_m2, dCR_m2, dCS_m2, dIMR_m2, dISR_m2, dIMS_m2, dISS_m2, dRR_m2, dRS_m2, dDR_m2, dDS_m2,
                     dU_f2, dCR_f2, dCS_f2, dIMR_f2, dISR_f2, dIMS_f2, dISS_f2, dRR_f2, dRS_f2, dDR_f2, dDS_f2, dN_to2))
  #results3 <- list(c(dCR_m2+dIMR_m2+dISR_m2+dCR_f2+dIMR_f2+dISR_f2))
  return(results2)
  })
}
# Calculate the min and max for each parameter
param_ranges2 <- sapply(parameters2, function(x) c(x * 0.9, x * 1.1))
# Convert the param_ranges to a matrix for easier handling in LHS
param_ranges_matrix2 <- matrix(unlist(param_ranges2), nrow = length(parameters2), byrow = TRUE, dimnames = list(names(parameters2), c("min", "max")))
params_range2 <- data.frame(t(param_ranges2)) 
# Define the parameter names and their respective ranges for sensitivity analysis
#ParameterLabelNames####
param_names2 <- c(
  "δ1",  
  "δ2",  
  "Discharge U[f]",
  "Discharge U[m]",
  "Discharge CR[f]",
  "Discharge CR[m]",
  "Discharge CS[f]",
  "Discharge CS[m]",
  "μ0",
  "μ1",
  "μ2",
  "μ3",
  "μ4",
  "μ5",
  "μ6",
  "Ψ[m]",
  "Ψ[f]",
  "Ψt[m]",
  "Ψt[f]",
  "c",
  "β1[m]",
  "β2[m]",
  "β1[f]",
  "β2[f]",
  "γ1",
  "γ2",
  "γ3",
  "γ4",
  "ω1d[m]",
  "ω1r[m]",
  "ω1d[f]",
  "ω1r[f]",
  "ω2d[m]",
  "ω2r[m]",
  "ω2d[f]",
  "ω2r[f]",
  "ω3d[m]",
  "ω3r[m]",
  "ω3d[f]",
  "ω3r[f]",
  "ω4d[m]",
  "ω4r[m]",
  "ω4d[f]",
  "ω4r[f]",
  "α1[m]",
  "α2[m]",
  "α1[f]",
  "α2[f]",
  "ε1",
  "ε2",
  "ζ3[m]",
  "ζ3[f]",
  "ζ1[m]",
  "ζ1[f]",
  "ζ2[m]",
  "ζ2[f]",
  "ζ4[m]",
  "ζ4[f]",
  "ν1[m]",
  "ν1[f]",
  "ν2[m]",
  "ν2[f]",
  "ν3[m]",
  "ν3[f]",
  "ν4[m]",
  "ν4[f]",
  "b",
  "φ[m]",
  "φ[f]",
  "π",
  "τ",
  "caIha", 
  "test"
)

#####
# Number of samples
n_samples <- 100
# Generate LHS samples
lhs_samples2 <- randomLHS(n = n_samples, k = length(parameters2))
# Scale the LHS samples to your parameters' ranges
scaled_samples2 <- apply(lhs_samples2, MARGIN = 2, FUN = function(x, range) {
  range[1] + (range[2] - range[1]) * x
}, range = t(params_range2))

simulation_results <- vector("list", n_samples)
times <- seq(from = 0, to = 365, by = 1)  # Adjust time range as needed

for (i in 1:n_samples) {
  parms_i <- setNames(as.list(scaled_samples2[i, ]), names(parameters2))
  parms_i <- unlist(list((parms_i)))
  simulation_results[[i]] <- ode(y = state2, times = times, func = ARB_model_2, parms = parms_i, method = "rk4")
}

#create outputs from above.
final_values_X1 <- numeric(n_samples)  # Replace 'n_samples' with the actual number of samples
final_values_X2 <- numeric(n_samples)  # Replace 'n_samples' with the actual number of samples
final_values_X3 <- numeric(n_samples)  # Replace 'n_samples' with the actual number of samples

for (i in 1:n_samples) {
  # Assuming 'X' is the name of your state variable and it's stored in a column named 'X'
  # Replace 'X' with the actual name or index of your variable
  final_values_X1[i] <- simulation_results[[i]][nrow(simulation_results[[i]]), "CR_m2"]+simulation_results[[i]][nrow(simulation_results[[i]]), "CR_f2"]+ simulation_results[[i]][nrow(simulation_results[[i]]), "IMR_m2"]+simulation_results[[i]][nrow(simulation_results[[i]]), "IMR_f2"]+ simulation_results[[i]][nrow(simulation_results[[i]]), "ISR_m2"]+simulation_results[[i]][nrow(simulation_results[[i]]), "ISR_f2"]
  final_values_X2[i] <- simulation_results[[i]][nrow(simulation_results[[i]]), "IMR_m2"]+simulation_results[[i]][nrow(simulation_results[[i]]), "IMR_f2"]+ simulation_results[[i]][nrow(simulation_results[[i]]), "ISR_m2"]+simulation_results[[i]][nrow(simulation_results[[i]]), "ISR_f2"]
  final_values_X3[i] <- simulation_results[[i]][nrow(simulation_results[[i]]), "DR_m2"]+simulation_results[[i]][nrow(simulation_results[[i]]), "DR_f2"]
}
results_df_allR <- data.frame(X_final_value = final_values_X1)
results_df_allRinf <- data.frame(X_final_value = final_values_X2)
results_df_deathsR <- data.frame(X_final_value = final_values_X3)

lhs_ranks <- apply(scaled_samples2 , 2, rank)
names(lhs_ranks)<-colnames(param_names2)
output_ranks <- apply(results_df_allR, 2, rank)
output_ranks <- as.numeric(output_ranks)  # Ensuring it's a numeric vector

library(corrr)
# Compute PRCC between each parameter and each output
prcc_results <- cor(lhs_ranks, output_ranks)
print(prcc_results)
prcc_values <- numeric(ncol(lhs_ranks))
for (i in 1:ncol(lhs_ranks)) {
  prcc_values[i] <- cor(lhs_ranks[, i], output_ranks, method = "spearman")
}
# Initialize a vector to hold p-values
prcc_p_values <- numeric(ncol(lhs_ranks))
# Perform the tests
for (i in 1:ncol(lhs_ranks)) {
  test_result <- cor.test(lhs_ranks[, i], output_ranks, method = "spearman")
  prcc_p_values[i] <- test_result$p.value
}
#ParameterLabelNames####
param_names2 <- c(
  "δ1",  
  "δ2",  
  "Discharge U[w]",
  "Discharge U[m]",
  "Discharge CR[w]",
  "Discharge CR[m]",
  "Discharge CS[w]",
  "Discharge CS[m]",
  "μ0",
  "μ1",
  "μ2",
  "μ3",
  "μ4",
  "μ5",
  "μ6",
  "Ψ[m]",
  "Ψ[w]",
  "Ψt[m]",
  "Ψt[w]",
  "c",
  "β1[m]",
  "β2[m]",
  "β1[w]",
  "β2[w]",
  "γ1",
  "γ2",
  "γ3",
  "γ4",
  "ω1d[m]",
  "ω1r[m]",
  "ω1d[w]",
  "ω1r[w]",
  "ω2d[m]",
  "ω2r[m]",
  "ω2d[w]",
  "ω2r[w]",
  "ω3d[m]",
  "ω3r[m]",
  "ω3d[w]",
  "ω3r[w]",
  "ω4d[m]",
  "ω4r[m]",
  "ω4d[w]",
  "ω4r[w]",
  "α1[m]",
  "α2[m]",
  "α1[w]",
  "α2[w]",
  "ε1",
  "ε2",
  "ζ3[m]",
  "ζ3[w]",
  "ζ1[m]",
  "ζ1[w]",
  "ζ2[m]",
  "ζ2[w]",
  "ζ4[m]",
  "ζ4[w]",
  "ν1[m]",
  "ν1[w]",
  "ν2[m]",
  "ν2[w]",
  "ν3[m]",
  "ν3[w]",
  "ν4[m]",
  "ν4[w]",
  "b",
  "φ[m]",
  "φ[w]",
  "π",
  "τ",
  "caIha", 
  "test"
)
#####
# Create the data frame using the named vectors
prcc_results2 <- data.frame(
  Parameter = param_names2,
  PRCC = prcc_values,
  P_Value = prcc_p_values  # Include this line only if you have computed p-values
)
prcc_results2$P_Value <- ifelse(prcc_results2$P_Value < 0.001, "<0.001", sprintf("%.3f", prcc_results2$P_Value))

#Latin Hypercube Sampling/Partial Rank Correlation Coefficient (LHS/PRCC) 
## Sort prcc_results based on the absolute values of PRCC but keep the original sign
prcc_results_sorted2 <- prcc_results2[order(abs(prcc_results2$PRCC), decreasing = TRUE), ]
library(ggplot2)
library(scales)
tiff("GlobalsensA_CRE.tiff", width = 13, height = 10, units = 'in', res = 500)
ggplot(prcc_results2, aes(x = reorder(Parameter, PRCC), y = PRCC, fill = PRCC)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flips the axes for easier reading
  scale_fill_gradient2(low = "cornflowerblue", mid = "beige", high = "gold2", midpoint = 0) +  # Color gradient
  geom_text(aes(label = paste("p ", P_Value), y = ifelse(PRCC < 0, PRCC - max(abs(PRCC) / 30), PRCC + max(abs(PRCC) / 30))), 
            position = position_dodge(width = 0.9), 
            hjust = ifelse(prcc_results2$PRCC < 0, 1.1, -0.1),  # Adjust text position based on PRCC value
            size = 3, family = "Times New Roman") +  # Adjust text size as needed
  labs(title = "",
       x = "Parameters",
       y = "PRCC value") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(-0.6, 0.5, by = 0.1), labels = scales::label_number(auto = TRUE))+ 
  theme(
  legend.position = "none",  # Remove legend if unnecessary
  axis.text.y = element_text(size = 10.5, family = "Times New Roman"),  # Set font family for y-axis text
  axis.text.x = element_text(family = "Times New Roman"),  # Set font family for x-axis text
  axis.title.y = element_text(family = "Times New Roman", margin = margin(r = 20)),  # Set font family for y-axis title
  axis.title.x = element_text(family = "Times New Roman"),  # Set font family for x-axis title
  text = element_text(family = "Times New Roman")  # Set font family for all text
)
dev.off()  # Close the


# Filter the results for significant p-values
significant_prcc_results2 <- prcc_results2 %>%
  filter(P_Value <= 0.05)
# Now plot using the filtered results
library(ggplot2)
# Plot using the filtered results
library(ggplot2)
# Plot using the filtered results
glob_sens_cre_pval <- ggplot(significant_prcc_results2, aes(x = reorder(Parameter, PRCC), y = PRCC, fill = PRCC)) +
  geom_bar(stat = "identity", color = "black") +  # Add black border around bars
  coord_flip() +  # Flips the axes for easier reading
  scale_fill_gradient2(low = "cornflowerblue", mid = "beige", high = "gold2", midpoint = 0) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) +  # Add horizontal line at y=0 (which is x=0 after flipping)
  geom_text(aes(label = paste("p = ", format(P_Value, digits = 2)), 
                y = ifelse(PRCC < 0, PRCC - max(abs(PRCC) / 30), PRCC + max(abs(PRCC) / 30))), 
            position = position_dodge(width = 0.1), 
            hjust = ifelse(significant_prcc_results2$PRCC < 0, 1.1, -0.1), 
            size = 5, family = "Times New Roman") +
  labs(title = "",
       x = "Parameters",
       y = "PRCC value") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 21, family = "Times New Roman"),
        axis.text.x = element_text(size = 21, family = "Times New Roman"),
        axis.title.y = element_text(family = "Times New Roman", size = 24, margin = margin(r = 20)),
        axis.title.x = element_text(family = "Times New Roman", size = 24),
        text = element_text(family = "Times New Roman")) +
  scale_y_continuous(breaks = seq(-0.8, 0.8, by = 0.2), labels = scales::label_number(auto = TRUE), limits = c(-0.9, 0.9))  # Extend y-axis limits
# Save the plot
ggsave("glob_sens_cre_pval.tiff", glob_sens_cre_pval, width = 12, height = 7, dpi = 800)

######################################################################################################################################################
####################################################################################################
######MIXING BOTH FIGURES:

# Load the necessary library
library(patchwork)
# Create a combined plot
combined_plot <- (glob_sens_cre_pval + annotate("text", x = Inf, y = -0.9, label = "(A) CRE", hjust = 0, vjust = 1, size = 6, fontface = "bold")) /
  (glob_sens_mrsa_pval + annotate("text", x = Inf, y = -0.9, label = "(B) MRSA", hjust = 0, vjust = 1, size = 6, fontface = "bold"))
combined_plot
# Save the combined plot
ggsave("combined_glob_sens.png", combined_plot, width = 13, height = 12, dpi = 800)

