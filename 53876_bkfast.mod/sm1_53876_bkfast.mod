NEURON
{
  SUFFIX bkfast 
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

  an = 0.1671594443234383     (/mV) 
  bn = -3.9301688070139904     (1) 
  vhn = -25.39201846411914     (mV) 
  An = 2.4322013030195464     (/ms) 
  b1n = -0.11709113141741097     (/mV) 
  c1n = 0.0015783247072133735     (/mV2) 
  d1n = -6.567272357576768e-06     (/mV3) 
  b2n = -0.07208406391765224     (/mV) 
  c2n = -0.00108305821104552     (/mV2) 
  d2n = -5.79690921581139e-06     (/mV3) 
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
  g = gbar*n*n
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