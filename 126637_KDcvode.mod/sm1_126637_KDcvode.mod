NEURON
{
  SUFFIX KDcvode 
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

  ah = -0.14417525969199135     (/mV) 
  bh = 13.692625630608129     (1) 
  vhh = -63.255682232349685     (mV) 
  Ah = 245.6767411295866     (/ms) 
  b1h = -0.014071997467490356     (/mV) 
  c1h = 0.0008977666565060346     (/mV2) 
  d1h = -4.728076010019702e-06     (/mV3) 
  b2h = 0.014071678747625622     (/mV) 
  c2h = -6.703926816710965e-05     (/mV2) 
  d2h = 2.625929411297013e-08     (/mV3) 

  am = 0.13766319121760132     (/mV) 
  bm = -6.170316068250969     (1) 
  vhm = -31.274609241639034     (mV) 
  Am = 8.615055796593243     (/ms) 
  b1m = -0.0894191904051355     (/mV) 
  c1m = 0.0014018686246833744     (/mV2) 
  d1m = -5.975540806783112e-06     (/mV3) 
  b2m = -0.007919423765154513     (/mV) 
  c2m = 0.0005086667960007121     (/mV2) 
  d2m = -2.7467023486406575e-06     (/mV3) 
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