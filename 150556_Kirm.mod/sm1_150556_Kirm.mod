NEURON
{
  SUFFIX Kirm 
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

  am = -0.0999998945510481     (/mV) 
  bm = 9.999993136233154     (1) 
  vhm = -506.7747749919132     (mV) 
  Am = 0.00503579929805654     (/ms) 
  b1m = 0.013921584794593449     (/mV) 
  c1m = -4.60837720491894e-05     (/mV2) 
  d1m = 4.474050405036246e-08     (/mV3) 
  b2m = -6.291518195640376     (/mV) 
  c2m = 0.021394323450133302     (/mV2) 
  d2m = -1.819287804629594e-05     (/mV3) 
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