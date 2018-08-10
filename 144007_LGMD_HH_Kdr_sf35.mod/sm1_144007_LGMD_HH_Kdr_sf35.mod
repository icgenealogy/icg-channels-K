NEURON
{
  SUFFIX HH_Kdr35 
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

  an = 0.09017259637329582     (/mV) 
  bn = -3.17896432391868     (1) 
  vhn = -45.892311033924095     (mV) 
  An = 0.7707743998875555     (/ms) 
  b1n = 0.06361500955593258     (/mV) 
  c1n = 0.0006064873584959908     (/mV2) 
  d1n = 4.7486223486674864e-06     (/mV3) 
  b2n = 0.038172420081031944     (/mV) 
  c2n = -0.00020995047320754848     (/mV2) 
  d2n = 4.815531218141136e-07     (/mV3) 
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
  g = gbar*n*n*n*n*n*n*n*n
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