NEURON
{
  SUFFIX ikdT 
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

  an = 0.09068612971092165     (/mV) 
  bn = -2.6114608696424333     (1) 
  vhn = -51.29893418383732     (mV) 
  An = 2.886790249536415     (/ms) 
  b1n = -0.04191407570642935     (/mV) 
  c1n = 0.0002739553551948287     (/mV2) 
  d1n = -7.539535344990324e-07     (/mV3) 
  b2n = -0.06990263673383365     (/mV) 
  c2n = -0.0011921993227104525     (/mV2) 
  d2n = -1.096148589889475e-05     (/mV3) 
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