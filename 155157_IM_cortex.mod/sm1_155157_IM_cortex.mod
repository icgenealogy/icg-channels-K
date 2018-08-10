NEURON
{
  SUFFIX im 
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

  am = 0.09996760613624721     (/mV) 
  bm = -3.49884126157024     (1) 
  vhm = -46.93922537213932     (mV) 
  Am = 506.4893847133652     (/ms) 
  b1m = 0.04999848792860153     (/mV) 
  c1m = -6.873214873875376e-08     (/mV2) 
  d1m = -1.17009361642987e-09     (/mV3) 
  b2m = 0.04999817256523821     (/mV) 
  c2m = 8.153415315503598e-08     (/mV2) 
  d2m = -1.119850993550486e-09     (/mV3) 
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