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

  am = 0.09984630685714171     (/mV) 
  bm = -2.9441063496114475     (1) 
  vhm = -10.374715451835783     (mV) 
  Am = 8.291511334071945     (/ms) 
  b1m = 0.18140937043119074     (/mV) 
  c1m = 0.0032628085973415897     (/mV2) 
  d1m = 1.919666185596494e-05     (/mV3) 
  b2m = 0.16351149807323906     (/mV) 
  c2m = -0.0024889829281940527     (/mV2) 
  d2m = 1.1918358169194043e-05     (/mV3) 
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
  g = gbar*m*m*m*m
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