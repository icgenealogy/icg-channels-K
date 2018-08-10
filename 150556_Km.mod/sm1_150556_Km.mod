NEURON
{
  SUFFIX Km 
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

  an = 0.057276757494883     (/mV) 
  bn = -1.3164033600502354     (1) 
  vhn = -49.57789609086683     (mV) 
  An = 2.3288508231040423     (/ms) 
  b1n = -0.03763697939256298     (/mV) 
  c1n = 0.00019480164102154073     (/mV2) 
  d1n = -4.1860349380397807e-07     (/mV3) 
  b2n = -0.03836162164586963     (/mV) 
  c2n = -0.0005000708771864327     (/mV2) 
  d2n = -3.66168465050199e-06     (/mV3) 
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