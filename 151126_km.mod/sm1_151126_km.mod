NEURON
{
  SUFFIX km 
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

  am = 0.25003713055341176     (/mV) 
  bm = -10.501740604625224     (1) 
  vhm = -74.72958970817409     (mV) 
  Am = 120.29996030968418     (/ms) 
  b1m = 0.01717110275642652     (/mV) 
  c1m = 0.0005365274903970418     (/mV2) 
  d1m = 1.5036203287220866e-05     (/mV3) 
  b2m = 0.015425053726481647     (/mV) 
  c2m = -0.00011413446612800195     (/mV2) 
  d2m = 2.795077016502037e-07     (/mV3) 
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