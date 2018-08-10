NEURON
{
  SUFFIX kdrgc 
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

  ah = -0.13424516877838186     (/mV) 
  bh = 7.032972970618415     (1) 
  vhh = -54.408076523685246     (mV) 
  Ah = 10424.486251046217     (/ms) 
  b1h = -0.09233802185662489     (/mV) 
  c1h = 0.0009587256553606102     (/mV2) 
  d1h = -3.1040155234834766e-06     (/mV3) 
  b2h = -0.08464829230068656     (/mV) 
  c2h = -0.0002008819027031214     (/mV2) 
  d2h = 4.287308942401115e-06     (/mV3) 

  am = 0.10638284513302973     (/mV) 
  bm = 0.021287809265012145     (1) 
  vhm = -29.91843539626418     (mV) 
  Am = 25.07586674799351     (/ms) 
  b1m = 0.08928276769968309     (/mV) 
  c1m = -1.8689402906978722e-05     (/mV2) 
  d1m = -3.1535673915155274e-07     (/mV3) 
  b2m = 0.022554726416870775     (/mV) 
  c2m = -4.210962340111969e-07     (/mV2) 
  d2m = -3.6309577109145673e-09     (/mV3) 
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
  g = gbar*h*m*m*m
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