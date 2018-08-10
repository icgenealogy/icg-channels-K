NEURON
{
  SUFFIX ichan2 
  USEION k READ ek WRITE ik 
  RANGE gbar, g, ik
  GLOBAL ek
}

UNITS
{
  (S) = (siemens)
  (mV) = (millivolt)
  (mA) = (milliamp)
}

PARAMETER
{
  gbar = 1 (S/cm2)

  ans = 0.12347619475896106     (/mV) 
  bns = -4.754174603420934     (1) 
  vhns = -45.325474447022984     (mV) 
  Ans = 0.5839592696636453     (/ms) 
  b1ns = 0.053320549702270195     (/mV) 
  c1ns = 0.0008200031204979145     (/mV2) 
  d1ns = 7.555821782871203e-06     (/mV3) 
  b2ns = 0.09288995599331366     (/mV) 
  c2ns = -0.0009188444084740952     (/mV2) 
  d2ns = 3.183772904054908e-06     (/mV3) 

  anf = 0.12349512899952127     (/mV) 
  bnf = -3.273305082351891     (1) 
  vhnf = -25.43918915703582     (mV) 
  Anf = 0.18398522406760032     (/ms) 
  b1nf = -0.10653614066809221     (/mV) 
  c1nf = 0.001976160275610243     (/mV2) 
  d1nf = -1.5247984159515117e-05     (/mV3) 
  b2nf = -0.017125092194956566     (/mV) 
  c2nf = 0.0001663843245489111     (/mV2) 
  d2nf = 1.5887675537941744e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  nsInf 
  nsTau 
  nfInf 
  nfTau 
}

STATE
{
  ns
  nf
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*ns*ns*ns*ns*nf*nf*nf*nf
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  ns' = (nsInf - ns) / nsTau 
  nf' = (nfInf - nf) / nfTau 
}

INITIAL
{
  rates(v)
  ns = nsInf 
  nf = nfInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nsInf = 1/(1 + exp(-ans*v + bns)) 
    nsTau = Ans / ( exp(-(b1ns*(v-vhns) + c1ns*(v-vhns)^2 + d1ns*(v-vhns)^3)) + exp((b2ns*(v-vhns) + c2ns*(v-vhns)^2 + d2ns*(v-vhns)^3)) ) 

    nfInf = 1/(1 + exp(-anf*v + bnf)) 
    nfTau = Anf / ( exp(-(b1nf*(v-vhnf) + c1nf*(v-vhnf)^2 + d1nf*(v-vhnf)^3)) + exp((b2nf*(v-vhnf) + c2nf*(v-vhnf)^2 + d2nf*(v-vhnf)^3)) ) 


  UNITSON
}