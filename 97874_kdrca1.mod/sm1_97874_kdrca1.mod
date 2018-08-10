NEURON
{
  SUFFIX kdr 
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

  an = 0.11222875641716407     (/mV) 
  bn = 1.459018644460495     (1) 
  vhn = 17.42806919482926     (mV) 
  An = 17.492308120931952     (/ms) 
  b1n = -0.10253784006498186     (/mV) 
  c1n = 0.0012981038995581206     (/mV2) 
  d1n = -2.4072965967235594e-06     (/mV3) 
  b2n = -0.01945530368058138     (/mV) 
  c2n = 0.000302687883212024     (/mV2) 
  d2n = 2.700035947244323e-06     (/mV3) 
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