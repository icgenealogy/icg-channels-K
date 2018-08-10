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

  ans = 0.12353245443998241     (/mV) 
  bns = -4.757007562037963     (1) 
  vhns = -43.73521434394381     (mV) 
  Ans = 0.5667011497154469     (/ms) 
  b1ns = 0.047115235827157174     (/mV) 
  c1ns = 0.0006339560120003629     (/mV2) 
  d1ns = 5.740937343371246e-06     (/mV3) 
  b2ns = 0.09544864112609074     (/mV) 
  c2ns = -0.0009764832990855176     (/mV2) 
  d2ns = 3.4734081639864975e-06     (/mV3) 

  anf = 0.12350333982435925     (/mV) 
  bnf = -3.2748044466084103     (1) 
  vhnf = -26.517597360131695     (mV) 
  Anf = 0.19377030557556038     (/ms) 
  b1nf = 0.022571817805983597     (/mV) 
  c1nf = -4.2812173968888174e-05     (/mV2) 
  d1nf = -7.194046886164697e-07     (/mV3) 
  b2nf = 0.106672798565709     (/mV) 
  c2nf = -0.001680547002817887     (/mV2) 
  d2nf = 8.516828015783e-06     (/mV3) 
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