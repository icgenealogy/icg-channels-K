NEURON
{
  SUFFIX iokdr 
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

  an = 0.05706747304794486     (/mV) 
  bn = -2.111028336968394     (1) 
  vhn = -61.31433817968842     (mV) 
  An = 1.1688513752996157     (/ms) 
  b1n = -0.04026424729087641     (/mV) 
  c1n = 0.00023993932913088872     (/mV2) 
  d1n = -6.035982768537479e-07     (/mV3) 
  b2n = -0.037502059899486566     (/mV) 
  c2n = -0.000478801714410455     (/mV2) 
  d2n = -3.3217076327334073e-06     (/mV3) 
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
  g = gbar*n*n*n*n
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