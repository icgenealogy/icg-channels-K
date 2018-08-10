NEURON
{
  SUFFIX km 
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

  an = 0.010820185943998209     (/mV) 
  bn = 1.4479014883221415     (1) 
  vhn = -34.83624324571596     (mV) 
  An = 33.843627312135354     (/ms) 
  b1n = 0.06378740137911912     (/mV) 
  c1n = 0.0007856882867620112     (/mV2) 
  d1n = 4.398042208048806e-06     (/mV3) 
  b2n = 0.04516663214336896     (/mV) 
  c2n = -0.0003054675574128596     (/mV2) 
  d2n = 8.669725517568978e-07     (/mV3) 
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