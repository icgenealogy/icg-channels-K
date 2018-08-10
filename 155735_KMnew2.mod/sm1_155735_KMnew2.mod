NEURON
{
  SUFFIX KMnew2 
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

  am = 0.10000084793378677     (/mV) 
  bm = -3.500101317787248     (1) 
  vhm = -47.178996894505865     (mV) 
  Am = 110.10176156386039     (/ms) 
  b1m = 0.05063672652063118     (/mV) 
  c1m = 1.314494620606467e-05     (/mV2) 
  d1m = 9.813702256059667e-08     (/mV3) 
  b2m = 0.0494344976169954     (/mV) 
  c2m = 9.453471074800518e-06     (/mV2) 
  d2m = -5.644436702899795e-08     (/mV3) 
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