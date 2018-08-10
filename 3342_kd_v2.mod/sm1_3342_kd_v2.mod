NEURON
{
  SUFFIX kd 
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

  an = 0.6760308436417274     (/mV) 
  bn = -31.171276737972565     (1) 
  vhn = -57.84156965958094     (mV) 
  An = 3.318844736579678     (/ms) 
  b1n = -0.028981245237038555     (/mV) 
  c1n = -0.0009390918654639228     (/mV2) 
  d1n = -7.262434750396567e-06     (/mV3) 
  b2n = 0.10582244692063245     (/mV) 
  c2n = -0.00019520827763873562     (/mV2) 
  d2n = -4.2358508352726904e-06     (/mV3) 
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