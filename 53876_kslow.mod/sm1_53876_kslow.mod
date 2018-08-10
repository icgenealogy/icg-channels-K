NEURON
{
  SUFFIX kslow 
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

  an = 0.19503278166308646     (/mV) 
  bn = -8.777138106284399     (1) 
  vhn = -39.057291544904906     (mV) 
  An = 428.97151415639786     (/ms) 
  b1n = -0.002407517414236337     (/mV) 
  c1n = -0.00045711193430009205     (/mV2) 
  d1n = 3.0087111759010966e-06     (/mV3) 
  b2n = 0.14900534683467767     (/mV) 
  c2n = -0.003478034944176158     (/mV2) 
  d2n = 1.7720924915134658e-05     (/mV3) 
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