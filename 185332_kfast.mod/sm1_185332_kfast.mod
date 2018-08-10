NEURON
{
  SUFFIX kfast 
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

  an = 0.10984330143912116     (/mV) 
  bn = 0.49359945169177544     (1) 
  vhn = 8.832668452608793     (mV) 
  An = 6.009821442315484     (/ms) 
  b1n = -0.08236250726284611     (/mV) 
  c1n = 0.0007693112617506492     (/mV2) 
  d1n = -2.9353757634900346e-06     (/mV3) 
  b2n = -0.02396344970534842     (/mV) 
  c2n = -0.0001153982570652951     (/mV2) 
  d2n = -2.550165987035497e-07     (/mV3) 
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