NEURON
{
  SUFFIX akfast 
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

  an = 0.21798650950529952     (/mV) 
  bn = -9.220756735668406     (1) 
  vhn = -36.96788472308139     (mV) 
  An = 4.044955934363971     (/ms) 
  b1n = -0.1560589329221162     (/mV) 
  c1n = 0.0038352310538892274     (/mV2) 
  d1n = -2.0040791257632154e-05     (/mV3) 
  b2n = 0.004898806095980935     (/mV) 
  c2n = 0.0004019088205330575     (/mV2) 
  d2n = -2.767920627615585e-06     (/mV3) 
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