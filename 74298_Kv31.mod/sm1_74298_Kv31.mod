NEURON
{
  SUFFIX Kv31 
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

  ap = 0.1111028289914337     (/mV) 
  bp = 0.009459063215699599     (1) 
  vhp = -16.43624042430112     (mV) 
  Ap = 16.808856967988294     (/ms) 
  b1p = -0.07603896676871215     (/mV) 
  c1p = 0.0005918619874085099     (/mV2) 
  d1p = -1.3902536974513495e-06     (/mV3) 
  b2p = -0.064650998983909     (/mV) 
  c2p = -0.0001313093816546784     (/mV2) 
  d2p = 2.6316324586540125e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  pInf 
  pTau 
}

STATE
{
  p
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*p
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  p' = (pInf - p) / pTau 
}

INITIAL
{
  rates(v)
  p = pInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    pInf = 1/(1 + exp(-ap*v + bp)) 
    pTau = Ap / ( exp(-(b1p*(v-vhp) + c1p*(v-vhp)^2 + d1p*(v-vhp)^3)) + exp((b2p*(v-vhp) + c2p*(v-vhp)^2 + d2p*(v-vhp)^3)) ) 


  UNITSON
}