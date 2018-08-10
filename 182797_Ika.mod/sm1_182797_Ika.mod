NEURON
{
  SUFFIX Ika 
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

  ah = -0.16666661359804444     (/mV) 
  bh = 13.083333029536076     (1) 
  vhh = 97.08292647060902     (mV) 
  Ah = 14.75296899623754     (/ms) 
  b1h = 0.004995791673843135     (/mV) 
  c1h = -5.093675150477357e-06     (/mV2) 
  d1h = 1.392298143849398e-07     (/mV3) 
  b2h = -0.004993570779244853     (/mV) 
  c2h = 5.075210293425605e-06     (/mV2) 
  d2h = -1.393177613751815e-07     (/mV3) 

  am = 0.049370616430092926     (/mV) 
  bm = 0.24100591933770146     (1) 
  vhm = -74.61998832140905     (mV) 
  Am = 0.22000000007231466     (/ms) 
  b1m = 6.236749821010841e-08     (/mV) 
  c1m = -3.523799547645604e-10     (/mV2) 
  d1m = 3.5082816464135165e-11     (/mV3) 
  b2m = 6.231729827335251e-08     (/mV) 
  c2m = -3.5017587322675286e-10     (/mV2) 
  d2m = 3.506587673046952e-11     (/mV3) 
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