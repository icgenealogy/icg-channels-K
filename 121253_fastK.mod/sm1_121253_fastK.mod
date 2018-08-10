NEURON
{
  SUFFIX fastK 
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

  an = 0.11536204578082107     (/mV) 
  bn = -7.364885082015897     (1) 
  vhn = -79.56226683597212     (mV) 
  An = 1.1553201973884015     (/ms) 
  b1n = -0.026075834557048164     (/mV) 
  c1n = 0.00011279449331845688     (/mV2) 
  d1n = -2.0287062675333388e-07     (/mV3) 
  b2n = -0.05221677097533658     (/mV) 
  c2n = -5.379232058815632e-05     (/mV2) 
  d2n = 1.7198154346639768e-06     (/mV3) 
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