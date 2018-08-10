NEURON
{
  SUFFIX Kdr_ChannelML 
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

  am = 0.09998596574962385     (/mV) 
  bm = 2.100017398857871     (1) 
  vhm = -49.02179470992445     (mV) 
  Am = 68.477909895746     (/ms) 
  b1m = -0.028168074674624913     (/mV) 
  c1m = 6.009432810152885e-06     (/mV2) 
  d1m = -1.771557011336049e-08     (/mV3) 
  b2m = -0.02670056047663398     (/mV) 
  c2m = 8.253520867402497e-06     (/mV2) 
  d2m = 2.3164533866453017e-08     (/mV3) 
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