NEURON
{
  SUFFIX Golgi_KM 
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

  an = 0.16667392093657826     (/mV) 
  bn = -5.833630896834246     (1) 
  vhn = 0.7778298179543348     (mV) 
  An = 48.434592410970325     (/ms) 
  b1n = -0.04751883029534987     (/mV) 
  c1n = 0.0006302190475145252     (/mV2) 
  d1n = -1.5755825426689038e-06     (/mV3) 
  b2n = -0.00660246077733082     (/mV) 
  c2n = 0.00023034955694828857     (/mV2) 
  d2n = -5.899635134853575e-07     (/mV3) 
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