NEURON
{
  SUFFIX GRC_KV 
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

  an = 0.05525123012910937     (/mV) 
  bn = -1.1452996119290122     (1) 
  vhn = -15.371012425134243     (mV) 
  An = 0.2901808556017063     (/ms) 
  b1n = -0.04653480362885211     (/mV) 
  c1n = 0.0005641250400951893     (/mV2) 
  d1n = -2.0459168581450346e-06     (/mV3) 
  b2n = -0.006032662851213341     (/mV) 
  c2n = 5.709509556987226e-05     (/mV2) 
  d2n = 2.186936350232098e-07     (/mV3) 
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