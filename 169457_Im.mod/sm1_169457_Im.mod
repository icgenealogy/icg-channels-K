NEURON
{
  SUFFIX Im 
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

  am = 0.20000755066221323     (/mV) 
  bm = -7.000343688665284     (1) 
  vhm = -35.020699990766815     (mV) 
  Am = 102.63399606429478     (/ms) 
  b1m = 0.10008322653555426     (/mV) 
  c1m = -2.0061801598537195e-06     (/mV2) 
  d1m = -1.0051064793274856e-07     (/mV3) 
  b2m = 0.0996746450690715     (/mV) 
  c2m = 1.63662314816597e-05     (/mV2) 
  d2m = -2.699301478206933e-07     (/mV3) 
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