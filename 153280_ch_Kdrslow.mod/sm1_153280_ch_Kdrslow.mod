NEURON
{
  SUFFIX ch_Kdrslow 
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

  an = 0.12327800725454417     (/mV) 
  bn = -4.74506387806256     (1) 
  vhn = -44.62073084336981     (mV) 
  An = 11.874216326284596     (/ms) 
  b1n = -0.09314336558361509     (/mV) 
  c1n = 0.0008718617039723995     (/mV2) 
  d1n = -2.9941798746483018e-06     (/mV3) 
  b2n = -0.04996235737011793     (/mV) 
  c2n = -0.0006856374802326956     (/mV2) 
  d2n = -6.045352655220736e-06     (/mV3) 
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