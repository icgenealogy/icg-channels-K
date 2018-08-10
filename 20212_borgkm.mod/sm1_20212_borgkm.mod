NEURON
{
  SUFFIX borgkm 
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

  am = 0.26188820075314434     (/mV) 
  bm = -14.404262101919594     (1) 
  vhm = -55.464029714377894     (mV) 
  Am = 457.34939544924924     (/ms) 
  b1m = 0.06481824752094927     (/mV) 
  c1m = 0.0006769206917693851     (/mV2) 
  d1m = 9.88747050926148e-06     (/mV3) 
  b2m = 0.22231550903298425     (/mV) 
  c2m = -0.002965690298107751     (/mV2) 
  d2m = 1.1305396316704671e-05     (/mV3) 
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