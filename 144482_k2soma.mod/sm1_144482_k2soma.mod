NEURON
{
  SUFFIX k2soma 
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

  ah = -0.0943267250374067     (/mV) 
  bh = 5.471023127353991     (1) 
  vhh = 6.4890315900549815     (mV) 
  Ah = 120.98911761126365     (/ms) 
  b1h = -0.0035872367112785367     (/mV) 
  c1h = 7.782857939433127e-06     (/mV2) 
  d1h = -8.915098763134379e-09     (/mV3) 
  b2h = -0.0035034372092122693     (/mV) 
  c2h = -5.013761163741049e-06     (/mV2) 
  d2h = 5.287301313911572e-10     (/mV3) 

  am = 0.05882447529466098     (/mV) 
  bm = -0.5882036842738817     (1) 
  vhm = -40.86633848688211     (mV) 
  Am = 77.16124151479154     (/ms) 
  b1m = -0.045292392281055305     (/mV) 
  c1m = 0.0002195112321659639     (/mV2) 
  d1m = -1.9957078017562646e-07     (/mV3) 
  b2m = -0.04759219403122255     (/mV) 
  c2m = 9.290300593580611e-07     (/mV2) 
  d2m = 2.60510438982881e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hInf 
  hTau 
  mInf 
  mTau 
}

STATE
{
  h
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*m
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}