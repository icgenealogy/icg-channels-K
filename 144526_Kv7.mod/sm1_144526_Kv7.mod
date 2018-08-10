NEURON
{
  SUFFIX Kv7 
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

  am = 0.08598549556144597     (/mV) 
  bm = -2.966001497884186     (1) 
  vhm = -34.55996886953395     (mV) 
  Am = 78.19696182168599     (/ms) 
  b1m = -0.05310346954933033     (/mV) 
  c1m = -5.126663489701884e-06     (/mV2) 
  d1m = 4.730273938509624e-08     (/mV3) 
  b2m = -0.03273185941631238     (/mV) 
  c2m = 1.8857726969968837e-07     (/mV2) 
  d2m = 6.634999231149526e-09     (/mV3) 
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