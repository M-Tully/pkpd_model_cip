@ -1,143 +0,0 @@
functions {
real two_cpt_oral_firstabs(real dose,
                             int t,
                             real tLag,
                             real tDose,
                             real CL,
                             real Vc,
                             real Q,
                             real Vp,
                             real ka) {
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
}
  
 
data {
  int<lower=1> nPt; // number of patients
  int<lower=1> N; // total # obs
  int<lower=0> t[N]; // time at each point
  int<lower=1> id[N]; // patient ID at each point
  int<lower=1> nObs[nPt]; // each patients # measurements  //
  int<lower=1> tmax[nPt]; // each patients max time
  real<lower=0> conc[N]; //simulated pk data for each measured time point
  int<lower=1> start[nPt];
  int<lower=1> end[nPt];
  real<lower=0>  tLag[nPt];
  real<lower=0>  tDose[nPt];
  real<lower=0>  dose[nPt];
  vector[5] Low; //lower parameter limits 
  vector[5] Upp; //upper parameter limits
}

parameters {
  real<lower=0> sig;
  vector[5] phiPop;
  
  // Inter-Individual variability
  cholesky_factor_corr[5] L;
  vector<lower = 0>[5] omega;
  matrix[5, nPt] etaStd;
}

transformed parameters{

 vector<lower=0>[sum(tmax)] conc_hat; // this is for every time point
 
  matrix[nPt, 5] phiInd;
  matrix[nPt, 5] thetaInd;//<lower=0>
  vector[5] thetaPop;//<lower=0>
  
  matrix[5, 5] Cor;
  matrix[5, 5] Cov;
  
  // Population average parameters 
   thetaPop = (Upp .* (1 - inv_logit(phiPop))) + (Low .* inv_logit(phiPop));
  
  // Individual parameters
  phiInd = (rep_matrix(phiPop, nPt) + diag_pre_multiply(omega, L * etaStd))';
  thetaInd = (rep_matrix(Upp, nPt)' .* (1 - inv_logit(phiInd))) + 
  (rep_matrix(Low, nPt)' .* inv_logit(phiInd));
  
  // Calculate correlation and covariance matrices
  Cor = L * L';
  Cov = quad_form_diag(Cor, omega);

for(a in 1:tmax[1]){
conc_hat[a] =  two_cpt_oral_firstabs(dose[1],
  a, // t
  tLag[1],
  tDose[1],
  thetaInd[1, 1], // Cl
  thetaInd[1, 2], // Vc
  thetaInd[1, 3], // Q
  thetaInd[1, 4], // Vp
  thetaInd[1, 5]  // ka
  );}

if (nPt > 1){
for (a in 2:nPt){ 
  for(c in (sum(tmax[1:(a-1)])+1):sum(tmax[1:a])){
conc_hat[c] =  two_cpt_oral_firstabs(dose[a],
 c - (sum(tmax[1:(a-1)])), // t
 tLag[a],
  tDose[a],
  thetaInd[a, 1], // Cl
  thetaInd[a, 2], // Vc
  thetaInd[a, 3], // Q
  thetaInd[a, 4], // Vp
  thetaInd[a, 5]  // ka
  );}}}
} 


model {
  
  // Population average
  phiPop ~ normal(0, 1);
  
  // Inter-individual variability 
  L ~ lkj_corr_cholesky(2);
  omega ~ normal(0, 1);
  to_vector(etaStd) ~ normal(0, 1);

for(c in 1:nObs[1]){
  conc[c] ~ normal(conc_hat[t[c]], sig);}
  
if (nPt > 1){
for (d in 2:nPt){
for (f in 1:nObs[d]){
conc[start[d]+(f-1)] ~ normal(conc_hat[sum(tmax[1:(d-1)])+t[sum(nObs[1:(d-1)])+f]], sig);}
}}
}

generated quantities{
real predict_conc[sum(tmax)];
  for (i in 1:sum(tmax)){
 predict_conc[i] = normal_rng(conc_hat[i], sig);}
}
