NEURON
{
  SUFFIX khva 
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

  an = 0.15986683781851088     (/mV) 
  bn = -3.102895309226248     (1) 
  vhn = -17.839621550773586     (mV) 
  An = 3.4233909926655524     (/ms) 
  b1n = -0.12165970171497206     (/mV) 
  c1n = 0.000492456504290907     (/mV2) 
  d1n = -5.47362147986132e-06     (/mV3) 
  b2n = -0.041204633928557555     (/mV) 
  c2n = 0.0001930616121853142     (/mV2) 
  d2n = 1.4717924063777469e-06     (/mV3) 
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