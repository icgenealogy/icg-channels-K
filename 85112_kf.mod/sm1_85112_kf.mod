NEURON
{
  SUFFIX kf 
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

  an = 0.11556605182150688     (/mV) 
  bn = -6.108906336788085     (1) 
  vhn = -18.188082847464106     (mV) 
  An = 4.467005959263763     (/ms) 
  b1n = -0.006175514938389175     (/mV) 
  c1n = -0.0001631094661147389     (/mV2) 
  d1n = 1.0512998232701227e-06     (/mV3) 
  b2n = 0.02549443520773997     (/mV) 
  c2n = -0.0003401888094243453     (/mV2) 
  d2n = 1.6297411905375324e-06     (/mV3) 
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