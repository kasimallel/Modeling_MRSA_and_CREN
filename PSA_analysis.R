####PSA ANALYSES:
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/AMR_BSI_BurdenChile/0_Article_AMR Transmission dynamics Chile & LMICs/0_BSIModelling/0_analysis/0_Figures_model")

library(parallel)
library(deSolve)

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
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*(((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1*(1-(sens_chrom_a*(1/(turn_chrom_a))*reduc_conpre_a)) )*(U_f1+U_m1))))/Nt1_spec) + b_p1*(r_v2)) 
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  FOC_cr_1 + ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
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
  FOC_cr_1 <-   FOC_cr_1 <- (((tau_p1*(1-c_p1)*(((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1*(1-(sens_chrom2_a*(1/(turn_chrom2_a))*reduc_conpre_a)) )*(U_f1+U_m1))))/Nt1_spec) + b_p1*(r_v2)) 
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  FOC_cr_1 + ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
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
  FOC_cr_1 <-   FOC_cr_1 <- (((tau_p1*(1-c_p1)*(((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1*(1-(sens_pcr_a*(1/(turn_pcr_a))*reduc_conpre_a)))*(U_f1+U_m1))))/Nt1_spec) + b_p1*(r_v2)) 
  FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  FOC_u_1 <-  FOC_cr_1 + ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
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
  tau_p1_rs <- (57.312+62.088)*tau_p1 /((1-0.6069)*(62.088)+0.6069*(57.312)) #(CR_m+CR_f1)*tau_p1 /((1-0.6069)*(CR_f1)+0.6069*(CR_m))
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions  0.3*Beta*women*uncolonised/(Nt) + 0.7*(1-Clevel)*Beta*(men)*uncolonised/Nt  Clevel= Coverage*efficacy*OR
  #FOC_cr_1 <- ((((1-0.6069)*tau_p1_rs*(1-c_p1)*((CR_f1+IMR_f1+ISR_f1))*(U_f1+U_m1)))/Nt1_spec)+(((0.6069*tau_p1_rs*(1-c_p1)*((CR_m1+IMR_m1+ISR_m1)*(1-reduc_conpre_a*or_HR_scenarMen_a))*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2)
  #FOC_cs_1 <- ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
  #FOC_u_1 <-  ((tau_p1*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)
  
  FOC_cr_1 <- (((tau_p1*(1-c_p1)*(1-reduc_conpre_a)*((CR_f1+CR_m1+IMR_m1+IMR_f1+ISR_m1+ISR_f1)*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2))
  FOC_cs_1 <- ((tau_p1*(1-reduc_conpre_a)*((CS_f1+CS_m1+IMS_m1+IMS_f1+ISS_m1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
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
  tau_p1_rs <- (57.312+62.088)*tau_p1 /((1-0.6069)*(62.088)+0.6069*(57.312)) #(CR_m+CR_f1)*tau_p1 /((1-0.6069)*(CR_f1)+0.6069*(CR_m))
  
  # DEFINITION OF THE FORCE OF INFECTION
  #FOC functions  0.3*Beta*women*uncolonised/(Nt) + 0.7*(1-Clevel)*Beta*(men)*uncolonised/Nt  Clevel= Coverage*efficacy*OR
  FOC_cr_1 <- ((((1-0.6069)*tau_p1_rs*(1-c_p1)*((CR_f1+IMR_f1+ISR_f1))*(U_f1+U_m1)))/Nt1_spec)+(((0.6069*tau_p1_rs*(1-c_p1)*(1-reduc_conpre_a*or_HR_scenarMen_a)*((CR_m1+IMR_m1+ISR_m1))*(U_f1+U_m1)))/Nt1_spec) + b_p1*(r_v2)
  FOC_cs_1 <- ((tau_p1*(1-reduc_conpre_a)*((CS_m1+IMS_m1+ISS_m1)*(U_f1+U_m1)))/Nt1_spec)  + ((tau_p1*((CS_m1+IMS_f1+ISS_f1)*(U_f1+U_m1)))/Nt1_spec)  
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
  FOC_cs_1 <- ((tau_p1*((CS_f1+IMS_f1+ISS_f1)*(1-reduc_conpre_a)*(U_f1+U_m1)))/Nt1_spec) + ((tau_p1*((CS_m1+IMS_m1+ISS_m1)*(U_f1+U_m1)))/Nt1_spec)
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













###################################
# Define the general structure for the loop and results storage
simulation_results <- list()  # Store results for each model
N_psa<-100
times <- seq(from=0, to=365, by = 1)  # Simulate over a year

#ParametersNEWMRSA####### 
# ----------------------------------------#
#Generate list of parameters with their distributions
params_d<- data.frame(
  delta1_p1 = rep(0.0016,N_psa),  
  delta2_p1 = rep(0.0016,N_psa),  
  Disch_U_f_p1 = rep(1/6,N_psa),
  Disch_U_m_p1 = rep(1/6,N_psa),
  Disch_CR_f_p1 = rep(1/6,N_psa),
  Disch_CR_m_p1 = rep(1/6,N_psa),
  Disch_CS_f_p1 = rep(1/6,N_psa),
  Disch_CS_m_p1 = rep(1/6,N_psa),
  mu0_p1 = rep(0.52,N_psa),  # % of women among U
  mu1_p1 = rep(0.52,N_psa),  # % of women among CR
  mu2_p1 = rep(0.52,N_psa),  # % of women among CS
  mu3_p1 = rep(0.3165,N_psa),  # % of women among IMR
  mu4_p1 = rep(0.3889,N_psa),  # % of women among ISR
  mu5_p1 = rep(0.3938,N_psa),  # % of women among IMS
  mu6_p1 = rep(0.4487,N_psa),  # % of women among IMS
  #Exposure to anbiotics; treatment amongst susceptible populations
  psi_m_p1 = rep(0.1474,N_psa), # % of  individuals exposed  to vancomycin/penicillin among males
  psi_w_p1 = rep(0.184,N_psa),  # % of  individuals exposed to vancomycin/penicillin among females
  #Percentage of people under treatment for MRSA decolonisation
  psi_mtr_p1=rep(0.109,N_psa), #IEAT is 1.35 times higher among CR; hence if treatment is psi_m_p1; psi_mtr_p1=psi_m_p1/1.35
  psi_wtr_p1=rep(0.136,N_psa),
  #Fitness cost. c reduces the transmission rate among resistant strains [1/unit time] [%].
  c_p1=rep(0.09,N_psa),
  #Progression to the development of infection from colonisation among CR and CS states. 
  beta1_m_p1 = rep((1/21)*0.26,N_psa),  # inverse of LOS plus progression from colonisation to infection among males CR
  beta2_m_p1 = rep((1/11)*0.099,N_psa),  # inverse of LOS plus progression from colonisation to infection among males CS
  beta1_f_p1 = rep((1/29)*0.26,N_psa), # inverse of LOS plus progression from colonisation to infection among females CR
  beta2_f_p1 = rep((1/14)*0.099,N_psa),  # inverse of LOS plus progression from colonisation to infection among females CS
  #Natural clearance of mild and severe infections among CR and CS states, respectively [1/unit time] [%].
  gamma1_p1 = rep(0.001,N_psa),  # Natural clearance among mild infections R
  gamma2_p1 = rep(0.001,N_psa),  # Natural clearance among severe infections R
  gamma3_p1 = rep(0.001,N_psa),  # Natural clearance among mild infections S
  gamma4_p1 = rep(0.001,N_psa), # Natural clearance among severe infections S
  #Mean time of infection considering length of hospital stays [1/length of hospital stay]. 
  omega1_d_m_p1 = rep((1/20),N_psa), #IMR patients who died, male
  omega1_r_m_p1 = rep((1/23),N_psa), #IMR patients who recovered, male
  omega1_d_f_p1 = rep((1/13),N_psa), #IMR patients who died, female
  omega1_r_f_p1 = rep((1/26),N_psa), #IMR patients who recovered, female
  omega2_d_m_p1 = rep((1/11),N_psa),  #ISR patients who died, male
  omega2_r_m_p1 = rep((1/18),N_psa),  #ISR patients who recovered, male
  omega2_d_f_p1 = rep((1/14),N_psa), #ISR patients who died, female 
  omega2_r_f_p1 = rep((1/19),N_psa), #ISR patients who recovered, female
  omega3_d_m_p1 = rep((1/12),N_psa), #IMS patients who died, male
  omega3_r_m_p1 = rep((1/11),N_psa), #IMS patients who recovered, male
  omega3_d_f_p1 = rep((1/20),N_psa), #IMS patients who died, female
  omega3_r_f_p1 = rep((1/16),N_psa), #IMS patients who recovered, female
  omega4_d_m_p1 = rep((1/11),N_psa), #ISS patients who died, male
  omega4_r_m_p1 = rep((1/19),N_psa), #ISS patients who recovered, male
  omega4_d_f_p1 = rep((1/14),N_psa), #ISS patients who died, female
  omega4_r_f_p1 = rep((1/17),N_psa), #ISS patients who recovered, female
  #Percentage of inpatients with CR or CS, respectively, progressing to severe infection in intensive care units [1/unit time] [%].
  alpha1_f_p1 = rep(0.4196,N_psa), #% patients with CR progressing to severe infection, males
  alpha2_f_p1 = rep(0.3438,N_psa), #% patients with CS progressing to severe infection, males
  alpha1_m_p1 = rep(0.3517,N_psa), #% patients with CR progressing to severe infection, females
  alpha2_m_p1 = rep(0.2826,N_psa), #% patients with CS progressing to severe infection, females
  #Progression from mild to severe infection from IMR and IMS, respectively [1/unit time] [%].
  epsilon1_p1 = rep(0.01,N_psa), #progression from IMR to ISR
  epsilon2_p1 = rep(0.01,N_psa), #progression from IMS to ISS
  #Mortality rates from infection. ζ1 and ζ2 are mortality rates from mild and severe resistant infections, respectively. ζ3 and ζ4 are from mild and severe susceptible infections, respectively [1/unit time] [%].
  zeta3_m_p1 = rep(0.231,N_psa),  # mortality rate from IMS  , male
  zeta3_f_p1 = rep(0.231*2.07,N_psa),  # mortality rate from IMS  , female
  zeta1_m_p1 = rep(0.231*1.01,N_psa),  # mortality rate from IMR, male
  zeta1_f_p1 = rep(0.231*1.22,N_psa),  # mortality rate from IMR, female
  zeta2_m_p1 = rep(0.231*2.32,N_psa),  # mortality rate from ISR, male
  zeta2_f_p1 = rep(0.231*2.25,N_psa),  # mortality rate from ISR, female
  zeta4_m_p1 = rep(0.231*1.10,N_psa),  # mortality rate from ISS , male
  zeta4_f_p1 = rep(0.231*2.34,N_psa),  # mortality rate from ISS , female
  #Recovery rates from infection, including IMR, ISR, IMS and ISS due to treatment received [1/unit time] [%].
  nu1_m_p1 = rep((1-zeta1_m_p1),N_psa), #recovery rates from IMR, males
  nu1_f_p1 = rep((1-zeta1_f_p1),N_psa), #recovery rates from IMR, females
  nu2_m_p1 = rep((1-zeta2_m_p1),N_psa), #recovery rates from ISR, males
  nu2_f_p1 = rep((1-zeta2_f_p1),N_psa), #recovery rates from ISR, females
  nu3_m_p1 = rep((1-zeta3_m_p1),N_psa), #recovery rates from IMS, males
  nu3_f_p1 = rep((1-zeta3_f_p1),N_psa), #recovery rates from IMS, females
  nu4_m_p1 =  rep((1-zeta4_m_p1),N_psa), #recovery rates from ISS, males
  nu4_f_p1 =  rep((1-zeta4_f_p1),N_psa),#recovery rates from ISS, females
  #Constant background rate that captures transmission from non-human sources, horizontal transmission, or de novo emergence [1/unit time] [number].
  b_p1 = rep(0.01,N_psa),
  # Percentage of people with resistant infections receiving inappropriate empirical antibiotic treatment [1/unit time] [%].
  phi_m_p1 = rep(0.459,N_psa),  # 
  phi_f_p1 = rep(0.413,N_psa),  # 
  #Factor of burden associated with inappropriate empirical antibiotic treatment and increased ICU admission among resistant infections [1/unit time] [%].
  pi_p1= rep(1.35,N_psa),
  #Transmission parameter {update this correspondingly after calibrating it with real data}
  #tau_p1= 0.03461113, #estimated 
  tau_p1=rep(0.2229124,N_psa),
  #community-acquired infection upon hospital admission rate
  caIha_p1= rep(0.001,N_psa),
  #percentage of people tested
  test_p1=rbeta(N_psa,0.20*100,100-0.2*100), 
  HR_perc1=rbeta(N_psa,0.2*100, 100-0.2*100),
  
  or_HR_scenar1_a=rgamma(N_psa, shape=25, scale=1.04/25),
  or_HR_scenarMen_a=rgamma(N_psa, shape=25, scale=1.37/25),
  ### ### ### ### ###
  #sensitivity chrom_1
  sens_chrom_a=rbeta(N_psa,0.826*100,100-0.826*100),
  #sensitivity chrom_1
  sens_chrom2_a=rbeta(N_psa,0.622*100,100-0.622*100), 
  #sensitivity chrom_1
  sens_pcr_a=rbeta(N_psa,0.881*100,100-0.881*100),
  #turnaround chrom_1
  turn_chrom_a=rgamma(N_psa, shape=25, scale=3/25),
  #turnaround chrom_1
  turn_chrom2_a=rgamma(N_psa, shape=25, scale=2/25),
  #turnaround pcr_1
  turn_pcr_a=rgamma(N_psa, shape=25, scale=1/25),
  #isolation contact precaution transmission reduction
  reduc_conpre_a=rbeta(N_psa,0.365*100, 100-0.365*100),
  #efficiency decolonisation
  eff_decol_a= rbeta(N_psa,0.53*100, 100-0.53*100),
  #effect on self-infection decolonisation
  eff_decol_selfi_a=rbeta(N_psa,0.32*100,100-0.32*100),
  #Turnaround decolonisation program in days
  turnaround_decol_a=rgamma(N_psa, shape=25, scale=5/25),
  #costs wards
  c_general_ward= rgamma(N_psa, shape=25, scale=50/25),
  c_intermediate_ward=rgamma(N_psa, shape=25, scale=92/25),
  c_icu_ward=rgamma(N_psa, shape=25, scale=218/25),
  c_decol_1pd=rgamma(N_psa, shape=25, scale=6.5/25),
  c_isolation=rgamma(N_psa, shape=25, scale=42.3/25),
  c_chrom=rgamma(N_psa, shape=25, scale=10.2/25),
  c_chrom2=rgamma(N_psa, shape=25, scale=13.6/25), 
  c_pcr=rgamma(N_psa, shape=25, scale=33/25),
  c_bc=rgamma(N_psa, shape=25, scale=16.9/25),
  #utilities
  u_healthy=rbeta(N_psa,0.92*100, 100-0.92*100),
  u_icu=rbeta(N_psa, 0.58*100,100-0.58*100),
  u_gw=rbeta(N_psa, 0.64*100,100-0.64*100),
  u_recovICU=rbeta(N_psa, 0.74*100, 100-0.74*100))

#####
model_list <- c("ARB_model_1ch_td_newadm","ARB_model_1ch2_td_newadm","ARB_model_1pcr_td_newadm","ARB_model_1ch_tiso_newadm","ARB_model_1ch2_tiso_newadm","ARB_model_1pcr_tiso_newadm",
                "ARB_model_1ch_td_newadmHR_m","ARB_model_1ch2_td_newadmHR_m","ARB_model_1_pcr_td_newadmHR_m","ARB_model_1ch_td_newadmHR_f","ARB_model_1ch2_td_newadmHR_f","ARB_model_1_pcr_td_newadmHR_f","ARB_model_1preE_newadm_all","ARB_model_1preE_newadm_m","ARB_model_1preE_newadm_f")  # Add last interventions with high-risk groups and pre-emptive isolation
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
#####
# ----------------------------------#
constant_cost_reduction<- results_icer[2,1]
constant_econ_factor<-  results_icer[3,1]
constant_pop_donoth<-results_icer[12, 1]

#ALL STRATEGIES+
results_qaly <- matrix(nrow = N_psa, ncol = length(model_list))
results_cost <- matrix(nrow = N_psa, ncol = length(model_list))


#Staphylococcus aureus bseline conditions[MRSA/MSSA]####### 
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
###CORES SIMULATIONS!
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterExport(cl, c("model_list", "N_psa", "params_d", "state1", "times", "constant_econ_factor", "constant_cost_reduction", "get", "ode","constant_pop_donoth" ,"ARB_model_1ch_td_newadm","ARB_model_1ch2_td_newadm","ARB_model_1pcr_td_newadm","ARB_model_1ch_tiso_newadm","ARB_model_1ch2_tiso_newadm","ARB_model_1pcr_tiso_newadm",
                    "ARB_model_1ch_td_newadmHR_m","ARB_model_1ch2_td_newadmHR_m","ARB_model_1_pcr_td_newadmHR_m","ARB_model_1ch_td_newadmHR_f","ARB_model_1ch2_td_newadmHR_f","ARB_model_1_pcr_td_newadmHR_f","ARB_model_1preE_newadm_all","ARB_model_1preE_newadm_m","ARB_model_1preE_newadm_f"))  # Export necessary objects to each worker
clusterEvalQ(cl, library(deSolve))  # Ensure each worker has required libraries

# Define function for parallel execution
run_model_parallel <- function(model_name) {
  model_results <- list(QALY = numeric(N_psa), Cost = numeric(N_psa))
  for (i in 1:N_psa) {
    model_func <- get(model_name, envir = .GlobalEnv)
    O_solution <- ode(y = state1, times = times, func = model_func, parms = params_d[i,], method = "rk4")
    recovered_pop <- O_solution[366, "RR_m1"] + O_solution[366, "RR_f1"] + O_solution[366, "RS_m1"] + O_solution[366, "RS_f1"]
    total_pop <- sum(O_solution[1:366, c("U_m1", "U_f1", "CR_m1", "CR_f1", "CS_m1", "CS_f1", "IMR_m1", "IMR_f1", "IMS_m1", "IMS_f1", "ISR_m1", "ISR_f1", "ISS_m1", "ISS_f1")])
    IMR_IMS_pop <- sum(O_solution[1:366, c("IMR_m1", "IMR_f1", "IMS_m1", "IMS_f1")])
    ISR_ISS_pop <- sum(O_solution[1:366, c("ISR_m1", "ISR_f1", "ISS_m1", "ISS_f1")])
    discharge_pop <- O_solution[366, "discharge"]
    # Calculations
    model_results$QALY[i] <- (((sum(O_solution[1:366, c("U_m1", "U_f1", "CR_m1", "CR_f1", "CS_m1", "CS_f1")])) ) * 0.92) +
      ((IMR_IMS_pop ) * 0.64) +
      ((ISR_ISS_pop) *  0.58) +(recovered_pop + discharge_pop) * 0.92 +
      ifelse((constant_pop_donoth - (total_pop + discharge_pop + recovered_pop)) > 0, (constant_pop_donoth - (total_pop + discharge_pop + recovered_pop)) * 0.92, 0) -
      constant_econ_factor
    model_results$Cost[i] <- O_solution[366, "cost"] - constant_cost_reduction
  }
  return(model_results)
}
# Run simulations in parallel
simulation_results <- parLapply(cl, model_list, run_model_parallel)
# Stop the cluster once done to free up resources
stopCluster(cl)
#12:50 Stopped at 18:50: Duration 7 hours approx.
simulation_results_mrsa <- simulation_results
######################################################################
######################################################################
#######GRAPHS!###############################################################
######################################################################
######################################################################

setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/AMR_BSI_BurdenChile/0_Article_AMR Transmission dynamics Chile & LMICs/0_BSIModelling/0_analysis/0_Figures_model/PSA")

model_listX <- c(
  "S1: T+D, agar1",
  "S2: T+D, agar2",
  "S3: T+D, pcr",
  "S4: T+I, agar1",
  "S5: T+I, agar2",
  "S6: T+I, pcr",
  "S7: T+D men, agar1",
  "S8: T+D men, agar2",
  "S9: T+D men, pcr",
  "S10: T+D women, agar1",
  "S11: T+D women, agar2",
  "S12: T+D women, pcr",
  "S13: Pre-emptive I, all",
  "S14: Pre-emptive I, men",
  "S15: Pre-emptive I, women")



names(simulation_results) <- model_listX 
strategy_names <- names(simulation_results)

simulation_resultsf <- simulation_results[!model_listX %in% c("S7: T+D men, agar1", "S8: T+D men, agar2", "S9: T+D men, pcr",
                                                                "S10: T+D women, agar1", "S11: T+D women, agar2", "S12: T+D women, pcr")]
names(simulation_resultsf) <- gsub("S13", "S7", names(simulation_resultsf))
names(simulation_resultsf) <- gsub("S14", "S8", names(simulation_resultsf))
names(simulation_resultsf) <- gsub("S15", "S9", names(simulation_resultsf))
strategy_names <- names(simulation_resultsf)

plot_listKI <- list()
for (strategy in strategy_names) {
  # Extract the strategy-specific results
  data1 <- data.frame(QALY = simulation_resultsf[[strategy]]$QALY,
                      Cost = simulation_resultsf[[strategy]]$Cost)
  # Create the plot
  p <- ggplot(data1, aes(x = QALY, y = Cost/10000)) +
    geom_point(shape = 21, color = "black", fill = "gold", size = 6) +  # Circle with border and fill
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 26, vjust = 0.5, hjust = 0.5),
          axis.text.y = element_text(size=26, angle = 0, vjust = 0.5, hjust = 1),
          plot.title = element_text(size = 37, face = "bold"),
          plot.subtitle = element_text(size = 37),
          axis.title.x = element_text(size = 26, vjust = -0.2),
          axis.title.y = element_text(size = 26, vjust = 1.5)) +
    labs(title = paste(""),
         subtitle = strategy,
         x = "Incremental QALYs",
         y = "Incremental Costs ($ in 10,000s)") +
    theme_lancet() +
    scale_y_continuous(labels = label_comma()) +
    scale_y_continuous(breaks = seq(-1000, 1000, by = 200)) +
    #scale_x_continuous(breaks = seq(0, 400, by = 50)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)
  plot_listKI[[strategy]] <- p
  # Save the plot with a strategy-specific filename
  file_name <- paste("psa_", gsub(" ", "_", tolower(strategy)), ".tiff", sep = "")
  ggsave(file_name, plot = p, width = 10, height = 12, dpi = 800)
}

combined_plotXOX <- wrap_plots(plot_listKI, ncol = 2)



updated_plot_list <- lapply(plot_listKI, function(p) {
  p + theme(
    axis.text.x = element_text(angle = 18, vjust = 0.0, hjust = 0.0, size = 14),  # Adjust text size as needed
    axis.text.y = element_text(size = 15, angle = 0, vjust = 0.5, hjust = 1),    # Adjust text size as needed
    plot.title = element_text(size = 20, face = "bold"),                         # Adjust title text size
    plot.subtitle = element_text(size = 18),                                     # Adjust subtitle text size
    axis.title.x = element_text(size = 17, vjust = -0.2),                        # Adjust x-axis title text size
    axis.title.y = element_text(size = 16, vjust = 1.5)                          # Adjust y-axis title text size
  )
})
# Combine updated plots
combined_plotXOX <- wrap_plots(updated_plot_list, ncol = 2)
  # Adjust ncol and nrow as needed to fit your layout
# Save the combined plot
ggsave("combined_qalys_all_mrsa.tiff", combined_plotXOX, width = 16, height = 21, dpi = 1000)
library(openxlsx)
write.xlsx(simulation_results, file = "simulations_mrsaCostS_qalys.xlsx", sheetName = "Sheet1", colNames = TRUE, rowNames = TRUE)




#WTP:
model_listX <- names(simulation_results)  # or a predefined list of model names if different
lancet_colors <- c("#E41A1C","#ffa819","#FC8D62","#984EA3","#d12b90","#F781BF", "#4DAF4A","#66C2A5","#A6D854","#2746ba","#8DA0CB","#377EB8", "#A65628","#d18b3b","#E5C494")

# Loop over indices to maintain alignment with lancet_colors
for (i in 1:length(model_listX)) {
  model <- model_listX[i]
  data1 <- data.frame(QALY = simulation_results[[model]]$QALY,
                      Cost = simulation_results[[model]]$Cost)
  data1$ICER <- data1$Cost / data1$QALY
  
  # Compute cost-effective percentages for various WTP thresholds
  wtp_thresholds <- seq(100, 16000, by = 100)
  cost_effective_percentages <- sapply(wtp_thresholds, function(wtp) {
    mean(data1$ICER <= wtp, na.rm = TRUE) * 100
  })
  plot_data <- data.frame(WTP = wtp_thresholds, Percentage = cost_effective_percentages)
  
  # Generate the plot
  p2 <- ggplot(plot_data, aes(x = WTP, y = Percentage)) +
    geom_line(color = "black", size = 5) +  # Underline for the contour
    geom_line(color = lancet_colors[i], size = 4) +  # Actual line color
    theme_minimal() +
    scale_x_continuous(labels = scales::comma, breaks = seq(0, 16000, by = 1000), limits = c(0, 16000)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
    labs(
      subtitle = paste(model),  # Added strategy info
      x = "Willingness-to-Pay Threshold ($)",
      y = "Percentage of cost-effective ICERs (ie, ICER<WTP)"
    ) +
    theme_lancet() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 19, face = "bold"),
      axis.title = element_text(size = 19),
      plot.subtitle = element_text(size = 19)
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)
  
  # Save the plot
  file_name <- paste("wtp_", gsub("[:\\[\\]\"+]", "", gsub(" ", "_", tolower(model))), ".tiff", sep = "")
  ggsave(file_name, plot = p2, width = 12, height = 10, dpi = 800)
}

############################################################
############################################################
###2: ALL GRAPHS WTP COMbINED IN ONE!
############################################################
############################################################

plot_data_list <- list()

# Iterate through each model strategy
for (i in 1:length(model_listX)) {
  model <- model_listX[i]
  # Extract corresponding data for current strategy
  data1 <- data.frame(QALY = simulation_results[[model]]$QALY,
                      Cost = simulation_results[[model]]$Cost)
  
  # Calculate ICER for each entry
  data1$ICER <- data1$Cost / data1$QALY
  
  # Initialize a list to collect percentage of cost-effective ICERs for different WTP thresholds
  percentages <- numeric(length(wtp_thresholds))
  
  # Calculate the percentage of cost-effective ICERs for each WTP threshold
  for (j in 1:length(wtp_thresholds)) {
    wtp = wtp_thresholds[j]
    percentages[j] <- mean(data1$ICER <= wtp, na.rm = TRUE) * 100  # Percentage of ICERs below current threshold
  }
  
  # Combine WTP thresholds with calculated percentages
  plot_data_list[[i]] <- data.frame(Strategy = rep(model, length(wtp_thresholds)),
                                    WTP = wtp_thresholds,
                                    Percentage = percentages)
}

# Combine all strategies' plot data into one dataframe for plotting
combined_plot_data <- do.call(rbind, plot_data_list)
combined_plot_data$Strategy <- factor(combined_plot_data$Strategy, levels = model_listX)


#FILTER INFORMATION TO ONLY 9 interventions:####
df_filteredx0x2 <- combined_plot_data[!grepl("S7|S8|S9|S10|S11|S12", combined_plot_data$Strategy), ]
df_filteredx0x2$Strategy <- gsub("S13", "S7", df_filteredx0x2$Strategy)
df_filteredx0x2$Strategy <- gsub("S14", "S8", df_filteredx0x2$Strategy)
df_filteredx0x2$Strategy <- gsub("S15", "S9", df_filteredx0x2$Strategy)

# Plot all strategies on one graph
wtp_1<-ggplot(df_filteredx0x2, aes(x = WTP, y = Percentage, group = Strategy, color = Strategy)) +
  geom_line(size = 2.3) +  # Draw lines
  scale_color_manual(values = lancet_colors) +  # Map lancet_colors to strategies
  theme_minimal() +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 16000, by = 1000), limits = c(0, 16000)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(30, 100), breaks = seq(30, 100, by = 10)) +
  labs(
    subtitle = "(B) MRSA model",
    x = "Willingness-to-pay threshold ($)",
    y = "Percentage of cost-effective ICERs (ie, ICER<WTP)",
    color = "Strategy"
  ) +
  theme_lancet()+
  theme(
    axis.text.x = element_text(size=15, angle = 45, hjust = 1),
    axis.text.y=element_text(size=17),
    plot.title = element_text(size = 19, face = "bold"),
    axis.title.y = element_text(size = 21),
    plot.subtitle = element_text(size=20),
    legend.text = element_text(size = 15),  # Increase legend text size
    legend.title = element_text(size = 18),  # Incr
    axis.title.x = element_text(size = 20, margin = margin(t = 20, r = 0, b = 0, l = 0)))+  # Increase top margin
  geom_hline(yintercept = 30, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)
# Save the combined plot
ggsave("combined_wtp_strategies_mrsa.tiff", plot=wtp_1, width = 13, height = 9, dpi = 1000)




library(openxlsx)
write.xlsx(combined_plot_data, file = "wtp_mrsa_data.xlsx", sheetName = "Sheet1", colNames = TRUE, rowNames = TRUE)







#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

####CRE:

#############################################################################
#############################################################################
#############################################################################
#############################################################################

####PSA ANALYSES:
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/AMR_BSI_BurdenChile/0_Article_AMR Transmission dynamics Chile & LMICs/0_BSIModelling/0_analysis/0_Figures_model/PSA2")

library(parallel)
library(deSolve)

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
  FOC_cs_2 <- ((tau_p2*(1-reduc_conpre_a)*((CS_f2+CS_m2+IMS_m2+IMS_f2+ISS_m2+ISS_f2)*(U_f2+U_m2)))/Nt2_spec2)  
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



###################################
# Define the general structure for the loop and results storage
simulation_results2 <- list()  # Store results for each model
N_psa<-100
times <- seq(from=0, to=365, by = 1)  # Simulate over a year
# ----------------------------------------#
#Generate list of parameters with their distributions
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
params_d <- data.frame(
  #clearance -natural- parameters
  delta1_p2 = rep(0.001, N_psa),  # Clearance value from CR state
  delta2_p2 = rep(0.001, N_psa),  # Clearance value from CS state
  #Dicharge rates from uncolonised and colonised
  Disch_U_f_p2 = rep(1/6, N_psa),
  Disch_U_m_p2 = rep(1/6, N_psa),
  Disch_CR_f_p2 = rep(1/6, N_psa),
  Disch_CR_m_p2 = rep(1/6, N_psa),
  Disch_CS_f_p2 = rep(1/6, N_psa),
  Disch_CS_m_p2 = rep(1/6, N_psa),
  #Proportion of women among specific populations (i.e., U, CR, CS, IMR, ISR, IMS, and ISS) [1/unit time] [%]
  mu0_p2 = rep(0.52, N_psa),  # % of women among U
  mu1_p2 = rep(0.52, N_psa),  # % of women among CR
  mu2_p2 = rep(0.52, N_psa),  # % of women among CS
  mu3_p2 = rep(0.3394, N_psa),  # % of women among IMR
  mu4_p2 = rep(0.3280, N_psa),  # % of women among ISR
  mu5_p2 = rep(0.50, N_psa),  # % of women among IMS
  mu6_p2 = rep(0.44, N_psa),  # % of women among IMS
  
  #Exposure to anbiotics 
  psi_m_p2 = rep(0.2225, N_psa), # % of  individuals exposed  to vancomycin/penicillin among males
  psi_w_p2 = rep(0.2026, N_psa),  # % of  individuals exposed to vancomycin/penicillin among males
  
  #Percentage of people under treatment for CRE decolonisation
  psi_mtr_p2=rep(0.2225/1.02, N_psa),
  psi_wtr_p2=rep(0.2026/1.02, N_psa),
  
  #Fitness cost. c reduces the transmission rate among resistant strains [1/unit time] [%].
  c_p2=rep((1-0.927), N_psa),
  
  #Progression to the development of infection from colonisation among CR and CS states. 
  beta1_m_p2 = rep((1/22)*0.213, N_psa), # inverse of LOS plus progression from colonisation to infection among males CR
  beta2_m_p2 = rep((1/20)*0.034, N_psa),  # inverse of LOS plus progression from colonisation to infection among males CS
  beta1_f_p2 = rep((1/27)*0.213, N_psa), # inverse of LOS plus progression from colonisation to infection among females CR
  beta2_f_p2 = rep((1/17)*0.034, N_psa),  # inverse of LOS plus progression from colonisation to infection among females CS
  
  #Natural clearance of mild and severe infections among CR and CS states, respectively [1/unit time] [%].
  gamma1_p2 = rep(0.001, N_psa),  # Natural clearance among mild infections R
  gamma2_p2 = rep(0.001, N_psa),  # Natural clearance among severe infections R
  gamma3_p2 = rep(0.001, N_psa),  # Natural clearance among mild infections S
  gamma4_p2 = rep(0.001, N_psa), # Natural clearance among severe infections S
  
  #Mean time of infection considering length of hospital stays [1/length of hospital stay]. 
  omega1_d_m_p2 =rep((1/21), N_psa), #IMR patients who died, male
  omega1_r_m_p2 =rep((1/26), N_psa), #IMR patients who recovered, male
  omega1_d_f_p2 =rep((1/31), N_psa), #IMR patients who died, female
  omega1_r_f_p2 =rep((1/30), N_psa), #IMR patients who recovered, female
  omega2_d_m_p2 =rep((1/7), N_psa), #ISR patients who died, male
  omega2_r_m_p2 =rep((1/20), N_psa),  #ISR patients who recovered, male
  omega2_d_f_p2 = rep((1/20), N_psa), #ISR patients who died, female 
  omega2_r_f_p2 =rep((1/23), N_psa), #ISR patients who recovered, female
  omega3_d_m_p2 =rep((1/12), N_psa), #IMS patients who died, male
  omega3_r_m_p2 =rep((1/20), N_psa), #IMS patients who recovered, male
  omega3_d_f_p2 =rep((1/10), N_psa), #IMS patients who died, female
  omega3_r_f_p2 =rep((1/18), N_psa), #IMS patients who recovered, female
  omega4_d_m_p2 =rep((1/11), N_psa), #ISS patients who died, male
  omega4_r_m_p2 =rep((1/14), N_psa), #ISS patients who recovered, male
  omega4_d_f_p2 =rep((1/9), N_psa), #ISS patients who died, female
  omega4_r_f_p2 =rep((1/15), N_psa), #ISS patients who recovered, female
  
  #Percentage of inpatients with CR or CS, respectively, progressing to severe infection in intensive care units [1/unit time] [%].
  alpha1_f_p2 = rep(0.4283, N_psa), #% patients with CR progressing to severe infection, males
  alpha2_f_p2 = rep(0.3548, N_psa), #% patients with CS progressing to severe infection, males
  alpha1_m_p2 = rep(0.4585, N_psa), #% patients with CR progressing to severe infection, females
  alpha2_m_p2 = rep(0.3832, N_psa), #% patients with CS progressing to severe infection, females
  
  #Progression from mild to severe infection from IMR and IMS, respectively [1/unit time] [%].
  epsilon1_p2 = rep(0.01, N_psa), #progression from IMR to ISR
  epsilon2_p2 = rep(0.01, N_psa), #progression from IMS to ISS
  
  #Mortality rates from infection. ζ1 and ζ2 are mortality rates from mild and severe resistant infections, respectively. ζ3 and ζ4 are from mild and severe susceptible infections, respectively [1/unit time] [%].
  zeta3_m_p2 = rep(0.228, N_psa),  # mortality rate from IMS  , male
  zeta3_f_p2 = rep(zeta3_m_p2*0.81, N_psa),  # mortality rate from IMS  , female
  zeta1_m_p2 = rep(zeta3_m_p2*1.80, N_psa),  # mortality rate from IMR, male
  zeta1_f_p2 = rep(zeta3_m_p2*0.55, N_psa),  # mortality rate from IMR, female
  zeta2_m_p2 = rep(zeta3_m_p2*1.30, N_psa),  # mortality rate from ISR, male
  zeta2_f_p2 = rep(zeta3_m_p2*2.40, N_psa),  # mortality rate from ISR, female
  zeta4_m_p2 = rep(zeta3_m_p2*1.62, N_psa),  # mortality rate from ISS , male
  zeta4_f_p2 = rep(zeta3_m_p2*2.23, N_psa),  # mortality rate from ISS , female
  
  #Recovery rates from infection, including IMR, ISR, IMS and ISS due to treatment received [1/unit time] [%].
  nu1_m_p2 = rep((1-zeta1_m_p2), N_psa), #recovery rates from IMR, males
  nu1_f_p2 = rep((1-zeta1_f_p2), N_psa), #recovery rates from IMR, females
  nu2_m_p2 = rep((1-zeta2_m_p2), N_psa), #recovery rates from ISR, males
  nu2_f_p2 = rep((1-zeta2_f_p2), N_psa), #recovery rates from ISR, females
  nu3_m_p2 = rep((1-zeta3_m_p2), N_psa), #recovery rates from IMS, males
  nu3_f_p2 = rep((1-zeta3_f_p2), N_psa), #recovery rates from IMS, females
  nu4_m_p2 =  rep((1-zeta4_m_p2), N_psa), #recovery rates from ISS, males
  nu4_f_p2 =  rep((1-zeta4_f_p2), N_psa),#recovery rates from ISS, females
  
  #Constant background rate that captures transmission from non-human sources, horizontal transmission, or de novo emergence [1/unit time] [number].
  b_p2 = rep(0.01, N_psa),
  
  # Percentage of people with resistant infections receiving inappropriate empirical antibiotic treatment [1/unit time] [%].
  phi_m_p2 = rep(0.3782, N_psa),  # Placeholder value, adjust as needed
  phi_f_p2 = rep(0.3846, N_psa),  # Placeholder value, adjust as needed
  
  #Factor of burden associated with inappropriate empirical antibiotic treatment and increased ICU admission among resistant infections [1/unit time] [%].
  pi_p2= rep(1.02, N_psa),
  
  #Transmission parameter {update this correspondingly after calibrating it with real data}
  tau_p2= rep(0.3986551, N_psa),
  
  #community-acquired infection upon hospital admission rate
  caIha_p2=rep(0.007, N_psa),
  
  #percentage of people tested
  test_p2=rbeta(N_psa, 0.20*100,100-0.2*100), 
  HR_perc2=rbeta(N_psa,0.2*100, 100-0.2*100),
  or_HR_scenar1_a=rgamma(N_psa, shape=25, scale=1.04/25),
  or_HR_scenarMen_a=rgamma(N_psa, shape=25, scale=2.27/25), 
  ### ### ### ### ###
  #sensitivity chrom_1
  sens_chrom_a=rbeta(N_psa,0.826*100,100-0.826*100),
  #sensitivity chrom_1
  sens_chrom2_a=rbeta(N_psa,0.90*100,100-0.9*100), 
  #sensitivity chrom_1
  sens_pcr_a=rbeta(N_psa,0.99*100, 100-0.99*100),
  #turnaround chrom_1
  turn_chrom_a=rgamma(N_psa, shape=25, scale=3/25),
  #turnaround chrom_1
  turn_chrom2_a=rgamma(N_psa, shape=25, scale=2/25),
  #turnaround pcr_1
  turn_pcr_a=rgamma(N_psa, shape=25, scale=1/25),
  #isolation contact precaution transmission reduction
  reduc_conpre_a=rbeta(N_psa,0.35*100,100-0.35*100),
  #efficiency decolonisation
  eff_decol_a= rbeta(N_psa,0.26*100,100-0.26*100), #0.146
  #effect on self-infection decolonisation
  eff_decol_selfi_a=rbeta(N_psa,0.041*100,100-0.041*100),
  #Turnaround decolonisation program in days
  turnaround_decol_a=rgamma(N_psa, shape=25, scale=7/25),
  ##
  #costs wards
  c_general_ward= rgamma(N_psa, shape=25, scale=50/25),
  c_intermediate_ward=rgamma(N_psa, shape=25, scale=92/25),
  c_icu_ward=rgamma(N_psa, shape=25, scale=218/25),
  c_decol_1pd=rgamma(N_psa, shape=25, scale=72.88/25),
  c_isolation=rgamma(N_psa, shape=25, scale=42.3/25),
  c_chrom= rgamma(N_psa, shape=25, scale=10.2/25),
  c_chrom2=rgamma(N_psa, shape=25, scale=13.6/25), 
  c_pcr=rgamma(N_psa, shape=25, scale=33/25),
  c_bc=rgamma(N_psa, shape=25, scale=16.9/25),
  #utilities
  u_healthy=rbeta(N_psa,0.92*100,100-0.92*100),
  u_icu=rbeta(N_psa,((0.92-0.34)*100),(100-(0.92-0.34)*100)),
  u_gw=rbeta(N_psa,(0.64*100),(100-100*0.64)),
  u_recovICU=rbeta(N_psa,0.74*100,100-0.74*100)
  
)

#####
model_list <- c("ARB_model_2ch_td_newadm","ARB_model_2ch2_td_newadm","ARB_model_2pcr_td_newadm","ARB_model_2ch_tiso_newadm","ARB_model_2ch2_tiso_newadm","ARB_model_2pcr_tiso_newadm",
                "ARB_model_2ch_td_newadmHR_m","ARB_model_2ch2_td_newadmHR_m","ARB_model_2_pcr_td_newadmHR_m","ARB_model_2ch_td_newadmHR_f","ARB_model_2ch2_td_newadmHR_f","ARB_model_2_pcr_td_newadmHR_f","ARB_model_2preE_newadm_all","ARB_model_2preE_newadm_m","ARB_model_2preE_newadm_f")  # Add last interventions with high-risk groups and pre-emptive isolation
#####
# ----------------------------------#
constant_cost_reduction<- results_icer[2,1]
constant_econ_factor<-  results_icer[3,1]
constant_pop_donoth<-results_icer[12, 1]
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

N_0m20 <-  U_m20 + CR_m20 + CS_m20 + IMR_m20 + ISR_m20 + IMS_m20 + ISS_m20 +  RR_m20 + RS_m20 + DR_m20 +DS_m20
N_0f20 <-  U_f20 + CR_f20 + CS_f20 + IMR_f20 + ISR_f20 + IMS_f20 + ISS_f20 +  RR_f20 + RS_f20 + DR_f20 +DS_f20
N_to0 <- N_0m20 + N_0f20

state2 <- c(U_m2 = U_m20, CR_m2=CR_m20, CS_m2= CS_m20, IMR_m2= IMR_m20, ISR_m2=ISR_m20, IMS_m2= IMS_m20, ISS_m2= ISS_m20, RR_m2= RR_m20, RS_m2=RS_m20, DR_m2= DR_m20, DS_m2=DS_m20,
            U_f2 = U_f20, CR_f2=CR_f20, CS_f2= CS_f20, IMR_f2= IMR_f20, ISR_f2=ISR_f20, IMS_f2= IMS_f20, ISS_f2= ISS_f20, RR_f2= RR_f20, RS_f2=RS_f20, DR_f2= DR_f20, DS_f2=DS_f20, N_to2=N_to0, utility=utility_to0, cost=cost_to0, new_admin=new_admin0, discharge=discharge0)
N_orig2<-N_0m20 + N_0f20
N_tdif <- N_orig2
Nt2_spec2<-1030


#ALL STRATEGIES+
results_qaly <- matrix(nrow = N_psa, ncol = length(model_list))
results_cost <- matrix(nrow = N_psa, ncol = length(model_list))
# Iterate over each model in the model list
#####
###CORES SIMULATIONS!
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterExport(cl, c("model_list", "N_psa","N_0m20","N_0f20", "params_d", "state2", "times", "constant_econ_factor", "constant_cost_reduction", "get", "ode","constant_pop_donoth" ,"ARB_model_2ch_td_newadm","ARB_model_2ch2_td_newadm","ARB_model_2pcr_td_newadm","ARB_model_2ch_tiso_newadm","ARB_model_2ch2_tiso_newadm","ARB_model_2pcr_tiso_newadm",
                    "ARB_model_2ch_td_newadmHR_m","ARB_model_2ch2_td_newadmHR_m","ARB_model_2_pcr_td_newadmHR_m","ARB_model_2ch_td_newadmHR_f","ARB_model_2ch2_td_newadmHR_f","ARB_model_2_pcr_td_newadmHR_f","ARB_model_2preE_newadm_all","ARB_model_2preE_newadm_m","ARB_model_2preE_newadm_f"))  # Export necessary objects to each worker
clusterEvalQ(cl, library(deSolve))  # Ensure each worker has required libraries
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

N_0m20 <-  U_m20 + CR_m20 + CS_m20 + IMR_m20 + ISR_m20 + IMS_m20 + ISS_m20 +  RR_m20 + RS_m20 + DR_m20 +DS_m20
N_0f20 <-  U_f20 + CR_f20 + CS_f20 + IMR_f20 + ISR_f20 + IMS_f20 + ISS_f20 +  RR_f20 + RS_f20 + DR_f20 +DS_f20
N_to0 <- N_0m20 + N_0f20

state2 <- c(U_m2 = U_m20, CR_m2=CR_m20, CS_m2= CS_m20, IMR_m2= IMR_m20, ISR_m2=ISR_m20, IMS_m2= IMS_m20, ISS_m2= ISS_m20, RR_m2= RR_m20, RS_m2=RS_m20, DR_m2= DR_m20, DS_m2=DS_m20,
            U_f2 = U_f20, CR_f2=CR_f20, CS_f2= CS_f20, IMR_f2= IMR_f20, ISR_f2=ISR_f20, IMS_f2= IMS_f20, ISS_f2= ISS_f20, RR_f2= RR_f20, RS_f2=RS_f20, DR_f2= DR_f20, DS_f2=DS_f20, N_to2=N_to0, utility=utility_to0, cost=cost_to0, new_admin=new_admin0, discharge=discharge0)
N_orig2<-N_0m20 + N_0f20
N_tdif <- N_orig2
Nt2_spec2<-1030
#####

#ALL STRATEGIES+
results_qaly <- matrix(nrow = N_psa, ncol = length(model_list))
results_cost <- matrix(nrow = N_psa, ncol = length(model_list))
# Iterate over each model in the model list

model_results <- list(QALY = numeric(N_psa), Cost = numeric(N_psa))
# Define function for parallel execution
run_model_parallel <- function(model_name) {
  model_results <- list(QALY = numeric(N_psa), Cost = numeric(N_psa))
  for (i in 1:N_psa) {
    #####
    model_func <- get(model_name, envir = .GlobalEnv)
    O_solution <- ode(y = state2, times = times, func = model_func, parms = params_d[i,], method = "rk4")
    recovered_pop <- O_solution[366, "RR_m2"] + O_solution[366, "RR_f2"] + O_solution[366, "RS_m2"] + O_solution[366, "RS_f2"]
    total_pop <- sum(O_solution[1:366, c("U_m2", "U_f2", "CR_m2", "CR_f2", "CS_m2", "CS_f2", "IMR_m2", "IMR_f2", "IMS_m2", "IMS_f2", "ISR_m2", "ISR_f2", "ISS_m2", "ISS_f2")])
    IMR_IMS_pop <- sum(O_solution[1:366, c("IMR_m2", "IMR_f2", "IMS_m2", "IMS_f2")])
    ISR_ISS_pop <- sum(O_solution[1:366, c("ISR_m2", "ISR_f2", "ISS_m2", "ISS_f2")])
    discharge_pop <- O_solution[366, "discharge"]
    # Calculations
    model_results$QALY[i] <- (((sum(O_solution[1:366, c("U_m2", "U_f2", "CR_m2", "CR_f2", "CS_m2", "CS_f2")])) ) * 0.92) +
      ((IMR_IMS_pop ) * 0.64) +
      ((ISR_ISS_pop) *  0.58) +(recovered_pop + discharge_pop) * 0.92 +
      ifelse((constant_pop_donoth - (total_pop + discharge_pop + recovered_pop)) > 0, (constant_pop_donoth - (total_pop + discharge_pop + recovered_pop)) * 0.92, 0) -
      constant_econ_factor
    model_results$Cost[i] <- O_solution[366, "cost"] - constant_cost_reduction
  }
  return(model_results)
}
# Run simulations in parallel
simulation_results2 <- parLapply(cl, model_list, run_model_parallel)
# Stop the cluster once done to free up resources
stopCluster(cl)
#12:50 Stopped at 18:50: Duration 7 hours approx.
simulation_resultsCRE<- simulation_results2



######################################################################
######################################################################
#######GRAPHS!###############################################################
######################################################################
######################################################################

setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/AMR_BSI_BurdenChile/0_Article_AMR Transmission dynamics Chile & LMICs/0_BSIModelling/0_analysis/0_Figures_model/PSA2")

model_listX <- c(
  "S1: T+D, agar1",
  "S2: T+D, agar2",
  "S3: T+D, pcr",
  "S4: T+I, agar1",
  "S5: T+I, agar2",
  "S6: T+I, pcr",
  "S7: T+D men, agar1",
  "S8: T+D men, agar2",
  "S9: T+D men, pcr",
  "S10: T+D women, agar1",
  "S11: T+D women, agar2",
  "S12: T+D women, pcr",
  "S13: Pre-emptive I, all",
  "S14: Pre-emptive I, men",
  "S15: Pre-emptive I, women")



names(simulation_results2) <- model_listX 
strategy_names <- names(simulation_results2)
simulation_results2f <- simulation_results2[!model_listX %in% c("S7: T+D men, agar1", "S8: T+D men, agar2", "S9: T+D men, pcr",
                                                        "S10: T+D women, agar1", "S11: T+D women, agar2", "S12: T+D women, pcr")]
names(simulation_results2f) <- gsub("S13", "S7", names(simulation_results2f))
names(simulation_results2f) <- gsub("S14", "S8", names(simulation_results2f))
names(simulation_results2f) <- gsub("S15", "S9", names(simulation_results2f))
strategy_names <- names(simulation_results2f)
plot_listKI <- list()
for (strategy in strategy_names) {
  # Extract the strategy-specific results
  data1 <- data.frame(QALY = simulation_results2f[[strategy]]$QALY,
                      Cost = simulation_results2f[[strategy]]$Cost)
  # Create the plot
  p <- ggplot(data1, aes(x = QALY, y = Cost/10000)) +
    geom_point(shape = 21, color = "black", fill = "gold", size = 6) +  # Circle with border and fill
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 26, vjust = 0.5, hjust = 0.5),
          axis.text.y = element_text(size=26, angle = 0, vjust = 0.5, hjust = 1),
          plot.title = element_text(size = 37, face = "bold"),
          plot.subtitle = element_text(size = 37),
          axis.title.x = element_text(size = 26, vjust = -0.2),
          axis.title.y = element_text(size = 26, vjust = 1.5)) +
    labs(title = paste(""),
         subtitle = strategy,
         x = "Incremental QALYs",
         y = "Incremental Costs ($ in 10,000s)") +
    theme_lancet() +
    scale_y_continuous(labels = label_comma()) +
    scale_y_continuous(breaks = seq(-1000, 1000, by = 200)) +
    #scale_x_continuous(breaks = seq(0, 400, by = 50)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)
  plot_listKI[[strategy]] <- p
  # Save the plot with a strategy-specific filename
  file_name <- paste("psa_", gsub(" ", "_", tolower(strategy)), ".tiff", sep = "")
  ggsave(file_name, plot = p, width = 10, height = 12, dpi = 800)
}

combined_plotXOX <- wrap_plots(plot_listKI, ncol = 2)



updated_plot_list <- lapply(plot_listKI, function(p) {
  p + theme(
    axis.text.x = element_text(angle = 18, vjust = 0.0, hjust = 0.0, size = 14),  # Adjust text size as needed
    axis.text.y = element_text(size = 15, angle = 0, vjust = 0.5, hjust = 1),    # Adjust text size as needed
    plot.title = element_text(size = 20, face = "bold"),                         # Adjust title text size
    plot.subtitle = element_text(size = 18),                                     # Adjust subtitle text size
    axis.title.x = element_text(size = 17, vjust = -0.2),                        # Adjust x-axis title text size
    axis.title.y = element_text(size = 16, vjust = 1.5)                          # Adjust y-axis title text size
  )
})
# Combine updated plots
combined_plotXOX <- wrap_plots(updated_plot_list, ncol = 2)
# Adjust ncol and nrow as needed to fit your layout
# Save the combined plot
ggsave("combined_qalys_all_cre.tiff", combined_plotXOX, width = 16, height = 21, dpi = 1000)
library(openxlsx)
write.xlsx(simulation_results2, file = "simulations_creCostS_qalys.xlsx", sheetName = "Sheet1", colNames = TRUE, rowNames = TRUE)




#WTP:
model_listX <- names(simulation_results2)  # or a predefined list of model names if different
#lancet_colors <- c("#E41A1C","#ffa819","#FC8D62","#984EA3","#E78AC3","#F781BF", "#4DAF4A","#66C2A5","#A6D854","#2746ba","#8DA0CB","#377EB8", "#A65628","#d18b3b","#E5C494")

# Loop over indices to maintain alignment with lancet_colors
for (i in 1:length(model_listX)) {
  model <- model_listX[i]
  data1 <- data.frame(QALY = simulation_results2[[model]]$QALY,
                      Cost = simulation_results2[[model]]$Cost)
  data1$ICER <- data1$Cost / data1$QALY
  
  # Compute cost-effective percentages for various WTP thresholds
  wtp_thresholds <- seq(100, 16000, by = 100)
  cost_effective_percentages <- sapply(wtp_thresholds, function(wtp) {
    mean(data1$ICER <= wtp, na.rm = TRUE) * 100
  })
  plot_data <- data.frame(WTP = wtp_thresholds, Percentage = cost_effective_percentages)
  
  # Generate the plot
  p2 <- ggplot(plot_data, aes(x = WTP, y = Percentage)) +
    geom_line(color = "black", size = 5) +  # Underline for the contour
    geom_line(color = lancet_colors[i], size = 4) +  # Actual line color
    theme_minimal() +
    scale_x_continuous(labels = scales::comma, breaks = seq(0, 16000, by = 1000), limits = c(0, 16000)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
    labs(
      subtitle = paste(model),  # Added strategy info
      x = "Willingness-to-Pay Threshold ($)",
      y = "Percentage of cost-effective ICERs (ie, ICER<WTP)"
    ) +
    theme_lancet() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 19, face = "bold"),
      axis.title = element_text(size = 19),
      plot.subtitle = element_text(size = 19)
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)
  
  # Save the plot
  file_name <- paste("wtp_", gsub("[:\\[\\]\"+]", "", gsub(" ", "_", tolower(model))), ".tiff", sep = "")
  ggsave(file_name, plot = p2, width = 12, height = 10, dpi = 800)
}

############################################################
############################################################
###2: ALL GRAPHS WTP COMbINED IN ONE!
############################################################
############################################################

plot_data_list <- list()

# Iterate through each model strategy
for (i in 1:length(model_listX)) {
  model <- model_listX[i]
  # Extract corresponding data for current strategy
  data1 <- data.frame(QALY = simulation_results2[[model]]$QALY,
                      Cost = simulation_results2[[model]]$Cost)
  
  # Calculate ICER for each entry
  data1$ICER <- data1$Cost / data1$QALY
  
  # Initialize a list to collect percentage of cost-effective ICERs for different WTP thresholds
  percentages <- numeric(length(wtp_thresholds))
  
  # Calculate the percentage of cost-effective ICERs for each WTP threshold
  for (j in 1:length(wtp_thresholds)) {
    wtp = wtp_thresholds[j]
    percentages[j] <- mean(data1$ICER <= wtp, na.rm = TRUE) * 100  # Percentage of ICERs below current threshold
  }
  
  # Combine WTP thresholds with calculated percentages
  plot_data_list[[i]] <- data.frame(Strategy = rep(model, length(wtp_thresholds)),
                                    WTP = wtp_thresholds,
                                    Percentage = percentages)
}

# Combine all strategies' plot data into one dataframe for plotting
combined_plot_data2 <- do.call(rbind, plot_data_list)
combined_plot_data2$Strategy <- factor(combined_plot_data2$Strategy, levels = model_listX)

# Plot all strategies on one graph
strategy_order <- c(
  "S1: T+D, agar1", "S2: T+D, agar2", "S3: T+D, pcr",
  "S4: T+I, agar1", "S5: T+I, agar2", "S6: T+I, pcr",
  "S7: T+D men, agar1", "S8: T+D men, agar2", "S9: T+D men, pcr",
  "S10: T+D women, agar1", "S11: T+D women, agar2", "S12: T+D women, pcr",
  "S13: Pre-emptive I, all", "S14: Pre-emptive I, men", "S15: Pre-emptive I, women"
)
combined_plot_data2 <- do.call(rbind, plot_data_list)
# Make sure 'Strategy' is a factor and order the levels as per strategy_order
combined_plot_data2$Strategy <- factor(combined_plot_data2$Strategy, levels = strategy_order)

#FILTER INFORMATION TO ONLY 9 interventions:####
df_filteredx0x <- combined_plot_data2[!grepl("S7|S8|S9|S10|S11|S12", combined_plot_data2$Strategy), ]
df_filteredx0x$Strategy <- gsub("S13", "S7", df_filteredx0x$Strategy)
df_filteredx0x$Strategy <- gsub("S14", "S8", df_filteredx0x$Strategy)
df_filteredx0x$Strategy <- gsub("S15", "S9", df_filteredx0x$Strategy)


# Plot all strategies on one graph
wtp_2 <- ggplot(df_filteredx0x, aes(x = WTP, y = Percentage, group = Strategy, color = Strategy)) +
  geom_line(size = 2.3) +  # Draw lines
  scale_color_manual(values = lancet_colors) +
  theme_minimal() +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 16000, by = 2000), limits = c(0, 16000)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(20, 100), breaks = seq(20, 100, by = 10)) +
  labs(
    subtitle = "(A) CRE model",
    x = "Willingness-to-pay threshold ($)",
    y = "Percentage of cost-effective ICERs",
    color = "Strategy"
  ) +
  theme_lancet() +
  theme(
    axis.text.x = element_text(size=15, angle = 45, hjust = 1),
    axis.text.y = element_text(size=17),
    plot.title = element_text(size = 19, face = "bold"),
    axis.title.y = element_text(size = 21),
    plot.subtitle = element_text(size=20),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 18),
    axis.title.x = element_text(size = 20, margin = margin(t = 20, r = 0, b = 0, l = 0))
  ) +
  geom_hline(yintercept = 20, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)

# Save the combined plot
ggsave("combined_wtp_strategies_cre.tiff", plot = wtp_2, width = 13, height = 9, dpi = 1000)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
# MERGING FIGURES CRE AND MRSA fOR WTP.
########################################################################################################################
########################################################################################################################
########################################################################################################################
lancet_colors <- c("#76b5c5","#f28e79","#a1c181","#f7bb5f","#4ebdb6","#ff9da7","#f4d4a4","#b699d7","#88d8b0")

library(cowplot)
# Combine wtp_1 and wtp_2 into one figure without their individual legends
wtp_1_updated <- wtp_1 + labs(y = "Percentage of cost-effective ICERs") + scale_x_continuous(breaks = seq(0, 16000, by = 2000), limits = c(0, 16000)) + theme_lancet()
wtp_2_updated <- wtp_2 + labs(y = "Percentage of cost-effective ICERs", x= "") +theme(axis.text.x = element_blank())


# Define a custom theme with Times New Roman font


combined_figure <- plot_grid(
  wtp_2_updated + theme_lancet()+geom_line(size = 2.8)+theme(legend.position = "none",  plot.subtitle = element_text(size=16), axis.text.x = element_blank()) ,  # Remove the legend from the first plot
  wtp_1_updated + theme_lancet()+geom_line(size = 2.8)+theme(legend.position = "none",  plot.subtitle = element_text(size=16), axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))  + scale_x_continuous(breaks = seq(0, 16000, by = 2000), limits = c(0, 16000)) ,  # Remove the legend from the second plot
  align = 'v',  # Align vertically
  #labels = c("(A)", "(B)"),  # Label your plots if needed
  ncol = 1  # Set number of columns
)

# Extract the legend from one of the plots (assuming both have the same legends)
legend <- get_legend(wtp_2)

# Combine the plots and the legend into a single figure, with the legend on the right
final_figure <- plot_grid(
  combined_figure,
  legend,
  ncol = 2,  # Place the legend as a second column
  rel_widths = c(30, 10)  # Adjust the relative widths between the plot and the legend
)

# Display the final combined figure
print(final_figure)
# Save the combined figure
ggsave("combined_wtp_strategies_MRSA_CRE.tiff", final_figure, width = 11.5, height = 13, dpi = 1000)

mean(simulation_results$`S2: T+D, agar2`$Cost)


