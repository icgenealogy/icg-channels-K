NEURON
{
  SUFFIX Kv1 
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

  an = 0.10999891299317557     (/mV) 
  bn = -4.949943318286745     (1) 
  vhn = -44.57891211702718     (mV) 
  An = 1.778699473352368     (/ms) 
  b1n = -0.03006924156405915     (/mV) 
  c1n = 5.1487804918708306e-06     (/mV2) 
  d1n = 3.4877297656412306e-08     (/mV3) 
  b2n = -0.07831939935451057     (/mV) 
  c2n = 0.00010756180546682958     (/mV2) 
  d2n = 2.892917805098672e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  nInf 
  nTau 
}

STATE
{
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}