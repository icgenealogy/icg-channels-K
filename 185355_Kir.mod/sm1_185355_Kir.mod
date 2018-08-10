NEURON
{
  SUFFIX kir 
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

  al = -0.09182682108894386     (/mV) 
  bl = 9.083536975974305     (1) 
  vhl = -89.48535880840024     (mV) 
  Al = 44.7289349497974     (/ms) 
  b1l = 0.01543304110826908     (/mV) 
  c1l = 3.604602771923279e-06     (/mV2) 
  d1l = 2.0323849241260124e-11     (/mV3) 
  b2l = 0.014349547368323522     (/mV) 
  c2l = 3.1113813396109286e-06     (/mV2) 
  d2l = -6.585761054076998e-09     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  lInf 
  lTau 
}

STATE
{
  l
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*l
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  l' = (lInf - l) / lTau 
}

INITIAL
{
  rates(v)
  l = lInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    lInf = 1/(1 + exp(-al*v + bl)) 
    lTau = Al / ( exp(-(b1l*(v-vhl) + c1l*(v-vhl)^2 + d1l*(v-vhl)^3)) + exp((b2l*(v-vhl) + c2l*(v-vhl)^2 + d2l*(v-vhl)^3)) ) 


  UNITSON
}