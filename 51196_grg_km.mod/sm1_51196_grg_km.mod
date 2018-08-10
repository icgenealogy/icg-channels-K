NEURON
{
  SUFFIX GrG_KM 
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

  an = 0.16667309536012734     (/mV) 
  bn = -5.000241384253188     (1) 
  vhn = -30.1533246172961     (mV) 
  An = 58.22704627872616     (/ms) 
  b1n = 0.050231468806516105     (/mV) 
  c1n = 4.089787468092437e-06     (/mV2) 
  d1n = 1.638117561943791e-08     (/mV3) 
  b2n = 0.0248103890737152     (/mV) 
  c2n = 2.1221820174796062e-06     (/mV2) 
  d2n = -9.362288280594827e-09     (/mV3) 
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