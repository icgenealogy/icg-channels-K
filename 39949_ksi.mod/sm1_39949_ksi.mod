NEURON
{
  SUFFIX ksi 
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

  an = 0.08474553032320394     (/mV) 
  bn = -1.144054807537135     (1) 
  vhn = -99.99988861145897     (mV) 
  An = 0.22268001412726024     (/ms) 
  b1n = -0.0032814282207581547     (/mV) 
  c1n = 2.9135845846881802e-05     (/mV2) 
  d1n = -7.847640735979825e-08     (/mV3) 
  b2n = 0.0032816167234399466     (/mV) 
  c2n = -2.9134161513436068e-05     (/mV2) 
  d2n = 7.847076693088014e-08     (/mV3) 
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