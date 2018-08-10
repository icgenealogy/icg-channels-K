NEURON
{
  SUFFIX Kdr 
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

  an = 0.06663263669145329     (/mV) 
  bn = -2.5316029284585104     (1) 
  vhn = -66.03522227428043     (mV) 
  An = 4.890416190086869     (/ms) 
  b1n = 0.025850284568278165     (/mV) 
  c1n = 5.1937991365665904e-05     (/mV2) 
  d1n = 3.367701532997099e-08     (/mV3) 
  b2n = 0.015224210214683412     (/mV) 
  c2n = 7.931548259487546e-05     (/mV2) 
  d2n = -2.28101928501713e-07     (/mV3) 
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