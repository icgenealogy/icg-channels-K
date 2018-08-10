NEURON
{
  SUFFIX KDR 
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

  an = 0.09719518906635281     (/mV) 
  bn = -1.6609222550081106     (1) 
  vhn = -35.00494672009407     (mV) 
  An = 10.20038660978946     (/ms) 
  b1n = 0.06571004928325934     (/mV) 
  c1n = 0.0009686626258843542     (/mV2) 
  d1n = 7.510929144783225e-06     (/mV3) 
  b2n = 0.05310903515442942     (/mV) 
  c2n = -0.0004208192331716756     (/mV2) 
  d2n = 1.3608385440654414e-06     (/mV3) 
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
  g = gbar*n
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