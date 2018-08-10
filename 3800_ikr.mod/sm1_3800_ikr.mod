NEURON
{
  SUFFIX IKr 
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

  am = 0.1538181299795938     (/mV) 
  bm = -2.1693863032970397     (1) 
  vhm = -17.96917599834666     (mV) 
  Am = 763.1376230554564     (/ms) 
  b1m = -0.0840991731129321     (/mV) 
  c1m = 0.0009132492733161869     (/mV2) 
  d1m = -3.771754434539875e-06     (/mV3) 
  b2m = -0.05189362044492394     (/mV) 
  c2m = -0.0006567712718443734     (/mV2) 
  d2m = -3.529718769288203e-06     (/mV3) 
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