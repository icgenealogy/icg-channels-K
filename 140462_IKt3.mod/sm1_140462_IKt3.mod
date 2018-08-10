NEURON
{
  SUFFIX Ikt3 
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

  an = 0.059847426085550194     (/mV) 
  bn = -1.5158800848929386     (1) 
  vhn = -32.834023079367206     (mV) 
  An = 1.3875083761686244     (/ms) 
  b1n = -0.03128851128372438     (/mV) 
  c1n = -0.0002010434297738282     (/mV2) 
  d1n = 2.090675190133158e-06     (/mV3) 
  b2n = -0.01082252848305164     (/mV) 
  c2n = 3.985595435580834e-05     (/mV2) 
  d2n = 1.086062562576792e-06     (/mV3) 
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