functions {

  // Initialise parasite age distribution  - truncated normal with only non-sequestered
  // this is because the initial time is inoculation, if not then should range over all ages
  vector PKPD_model_initialisation(int nStages, 
  real ipl, 
  real iplMu, 
  real iplSigma,
  real pmf,
  int tseq) {
    
    vector[nStages] x0; 
    vector[nStages] pAge;
    
     for (i in 1:(tseq)) { // 
      pAge[i] = (normal_lpdf(i | iplMu, iplSigma));
    }
    for ( i in (tseq+1):nStages){
      pAge[i] = 0;
    }
    
      x0[1:tseq] = ipl* softmax(pAge[1:tseq]);

    for ( i in (tseq+1):nStages){
      x0[i] = 0;
    }
    
    
    
    return x0; 
  }
  
  
  
 real two_cpt_oral_firstabs2(real dose,
                             int t,
                             real tLag,
                             real tDose,
                             real CL,
                             real Vc,
                             real Q,
                             real Vp,
                             real ka
                             ) {
    real beta;
    real alpha;
    real A;
    real B;
    real Ae;
    real Be;
    real Ce;
    real time;
    real conc;
  
  beta = (((CL/Vc)+(Q/Vc)+(Q/Vp)) - sqrt((((CL/Vc)+(Q/Vc)+(Q/Vp))^2)-4*(CL/Vc)*(Q/Vp)))/2;
  alpha = ((Q/Vp) * (CL/Vc)) / beta;
  
  A = (ka/Vc) * (Q/Vp - alpha)/((ka - alpha) * (beta - alpha));
  B = (ka/Vc) * (Q/Vp - beta)/((ka - beta) * (alpha - beta));
  Ae = A;
  Be = B;
  Ce = (Ae*(ka - alpha) + Be*(ka - beta)) / (ka);
  time = t - tLag - tDose;
  if (time < 0){
    conc = 0;
  } else {
    conc = 1000 * dose * (Ae * exp(-alpha*time) + Be * exp(-beta*time) + Ce * exp(-time) - (Ae+Be+Ce) * exp(-ka*time));
  if(conc < 0){ conc = 0;}
  }
  return conc;
}
  
  
  
  // model 
matrix parasite_mod(int nStages,
    real ipl, 
    real iplMu, 
    real iplSigma,
    real PMF,
    int tmax,
    real Emax,
    int nObs,
    real dose,
  real tLag,
  real tDose,
  real CL,
  real Vc,
  real Q,
  real Vp,
  real ka,
    real Ec50,
    real gamma,
    real[] kill_window,
    int tseq
    ){
      
    matrix[nStages,tmax] totCountAgeMatrix;
    vector[tmax]  NCount;
    vector[tmax]  cirCount;
    vector[tmax]  NCount_hat;
    vector[tmax]  cirCount_hat;
  
    vector[tmax] logNCount_hat;
    vector[tmax] logcirCount_hat;
    
    vector[tmax] E;
    matrix[nStages,tmax] survProb;
    
    matrix[2,tmax] result;
    
    vector[tmax] conc;
  
  for (x in 1:nStages){
  totCountAgeMatrix[x,1] =  PKPD_model_initialisation(nStages, 
   ipl, 
   iplMu, 
   iplSigma,
   PMF,
   tseq)[x];}

  NCount[1] = sum(totCountAgeMatrix[,1]);
  cirCount[1] = sum(totCountAgeMatrix[1:tseq,1]);
  

for (a in 1:tmax){
  conc[a] = two_cpt_oral_firstabs2(dose,
    a,
    tLag,
    tDose,
    CL,
    Vc,
    Q,
    Vp,
    ka);}
   for (i in 1:tmax){
    E[i] = (Emax * conc[i]^gamma) / (conc[i]^gamma + Ec50^gamma);
 
  for (a in 1:nStages){
   survProb[a,i] = 1 - E[i] * kill_window[a];}}
   
  for (j in 2:tmax){ 
  for (a in 2:nStages){
      totCountAgeMatrix[a,j] = totCountAgeMatrix[a-1,j-1] * survProb[a-1,j-1];}
      totCountAgeMatrix[1,j] = PMF * totCountAgeMatrix[nStages,j-1] * survProb[nStages,j-1]; 
   if (sum(totCountAgeMatrix[,j]) < 1) {
      for (h in 1:nStages){
      totCountAgeMatrix[h,j] = 0;} }
     NCount[j] = sum(totCountAgeMatrix[,j]);
     cirCount[j] = sum(totCountAgeMatrix[1:(tseq),j]); }
     
  for (p in 1:tmax){
     NCount_hat[p] = NCount[p];
     cirCount_hat[p] = cirCount[p];}
     
for (i in 1:tmax){ 
  
 if (NCount_hat[i] > 1){
  logNCount_hat[i] = log(NCount_hat[i]); }
  else {  logNCount_hat[i] = 0;}
  
 if (cirCount_hat[i] > 1){
  logcirCount_hat[i] = log(cirCount_hat[i]); }
  else {  logcirCount_hat[i] = 0;}

result[1,i] = logNCount_hat[i];
result[2,i] = logcirCount_hat[i];}

return result;
} 

}

data {
  int<lower=1> nStages; // i.e 48
  int<lower=1> nPt; // number of patients
  int<lower=1> N; // total # obs
  int<lower=0> t[N]; // time at each poin
  int<lower=1> nObs[nPt]; // each patients # measurements  //
  int<lower=1> tmax[nPt]; // each patients max time
  int<lower=1> id[sum(tmax)]; // patient ID at each point (long ver)
  int<lower=1> ID[N]; // patient ID at each obs point (short ver)
 // real<lower=0> conc[N]; //fitted pk data for each time point (made with input params)
  real<lower=0, upper=1> kill_window[nStages]; // binary vector
//  real<lower=0> NCount_measured_e[N];
  real<lower=0> cirCount_measured_e[N];
  int<lower=1> start[nPt];
  int<lower=1> end[nPt];
  real<lower=0>  tLag[nPt];
  real<lower=0>  tDose[nPt];
  real<lower=0>  dose[nPt];
  real<lower=0>  CL[nPt]; // will be estimates!
  real<lower=0>  Vc[nPt];
  real<lower=0>  Q[nPt];
  real<lower=0>  Vp[nPt];
  real<lower=0>  ka[nPt];
  vector[7] Low; //lower parameter limits 
  vector[7] Upp; //upper parameter limits
  int<lower=12> tseq;
}
parameters {

  real<lower=0,upper=5> se[nPt]; 
 
  vector[7] phiPop;
  
  // Inter-Individual variability
  cholesky_factor_corr[7] L;
  vector<lower = 0>[7] omega;
  matrix[7, nPt] etaStd;
}

transformed parameters{

  matrix[nPt, 7] phiInd;
  matrix[nPt, 7] thetaInd;//<lower=0>
  vector[7] thetaPop;//<lower=0>
  
  matrix[7, 7] Cor;
  matrix[7, 7] Cov;
  
   matrix<lower=0>[2,sum(tmax)] hat_matrix;
  
  // Population average parameters 
   thetaPop = (Upp .* (1 - inv_logit(phiPop))) + (Low .* inv_logit(phiPop));
  
  // Individual parameters
  phiInd = (rep_matrix(phiPop, nPt) + diag_pre_multiply(omega, L * etaStd))';
  thetaInd = (rep_matrix(Upp, nPt)' .* (1 - inv_logit(phiInd))) + 
  (rep_matrix(Low, nPt)' .* inv_logit(phiInd));
  
  // Calculate correlation and covariance matrices
  Cor = L * L';
  Cov = quad_form_diag(Cor, omega);


hat_matrix[,1:tmax[1]] = parasite_mod(nStages,
  thetaInd[1,1], //ipl[1], 
  thetaInd[1,2], // iplMu[1], 
  thetaInd[1,3], //iplSigma[1],
  thetaInd[1,4], //PMF[1],
  tmax[1],
  thetaInd[1,5], //Emax
  nObs[1],
  dose[1],
  tLag[1],
  tDose[1],
  CL[1],
  Vc[1],
  Q[1],
  Vp[1],
  ka[1],
  thetaInd[1,6], //Ec50
  thetaInd[1,7], //gamma
  kill_window,
  tseq);


if (nPt > 1){
for (a in 2:nPt){ 
hat_matrix[,(sum(tmax[1:((a-1))])+1):sum(tmax[1:a])] = parasite_mod(nStages,
  thetaInd[a,1], //ipl[1], 
  thetaInd[a,2], // iplMu[1], 
  thetaInd[a,3], //iplSigma[1],
  thetaInd[a,4], //PMF[1],
  tmax[a],
  thetaInd[a,5], //Emax
  nObs[a],
  dose[a],
  tLag[a],
  tDose[a],
  CL[a],
  Vc[a],
  Q[a],
  Vp[a],
  ka[a],
  thetaInd[a,6], //Ec50
  thetaInd[a,7], //gamma
  kill_window,
  tseq);
}}}


model {
  
  // Population average
  phiPop ~ normal(0, 1);
  
  // Inter-individual variability 
  L ~ lkj_corr_cholesky(2);
  omega ~ normal(0, 1);
  to_vector(etaStd) ~ normal(0, 1);


for (c in 1:N){

if (cirCount_measured_e[c] >= log(50000)) {
  // non=censored
cirCount_measured_e[c] ~ normal(hat_matrix[2,sum(tmax[1:(ID[c]-1)])+t[c]], se[ID[c]]);}
else{
  // censored
  target += normal_lcdf(cirCount_measured_e[c] | hat_matrix[2,sum(tmax[1:(ID[c]-1)])+t[c]], se[ID[c]]);}}
  
}

generated quantities{
real predict_cir[sum(tmax)];
real predict_tot[sum(tmax)];
  for (i in 1:(sum(tmax))){
   predict_tot[i] = normal_rng(hat_matrix[1,i], se[id[i]]);
  predict_cir[i] = normal_rng(hat_matrix[2,i], se[id[i]]);}
}
