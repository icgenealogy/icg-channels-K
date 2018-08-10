NEURON
{
  SUFFIX kv72wt73wt 
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

  am = 0.08583578009602673     (/mV) 
  bm = -2.6351569275191933     (1) 
  vhm = -8.772704839202817     (mV) 
  Am = 18.986231465565     (/ms) 
  b1m = 0.11234425134377485     (/mV) 
  c1m = 0.007570932992588774     (/mV2) 
  d1m = 6.831443075254374e-05     (/mV3) 
  b2m = -0.0011471792584379003     (/mV) 
  c2m = 0.0001685834744744023     (/mV2) 
  d2m = -1.090939659672431e-06     (/mV3) 
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