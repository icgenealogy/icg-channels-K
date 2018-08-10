NEURON
{
  SUFFIX kdrsup 
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

  am = 0.08264389775701347     (/mV) 
  bm = 0.2892668635911771     (1) 
  vhm = -48.59024551276862     (mV) 
  Am = 54.95422461092577     (/ms) 
  b1m = 0.026283077277156945     (/mV) 
  c1m = -1.5206091760122105e-05     (/mV2) 
  d1m = -7.858762201348081e-08     (/mV3) 
  b2m = 0.028414087498025986     (/mV) 
  c2m = -8.156484140663646e-06     (/mV2) 
  d2m = 2.4605996123570167e-08     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
}

STATE
{
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}