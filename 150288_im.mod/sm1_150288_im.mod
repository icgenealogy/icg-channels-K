NEURON
{
  SUFFIX im 
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

  an = 0.09667202654104037     (/mV) 
  bn = -5.0946628097001385     (1) 
  vhn = -53.02177941901189     (mV) 
  An = 62.40565959665734     (/ms) 
  b1n = -0.04283498829839039     (/mV) 
  c1n = -9.529962863905686e-06     (/mV2) 
  d1n = 5.345960717232205e-08     (/mV3) 
  b2n = -0.05405378990217072     (/mV) 
  c2n = -2.0149093883193045e-05     (/mV2) 
  d2n = -1.79960675590407e-07     (/mV3) 
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