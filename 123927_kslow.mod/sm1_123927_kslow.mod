NEURON
{
  SUFFIX kslow 
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

  ah = -0.09077475909704454     (/mV) 
  bh = 4.904797497115002     (1) 
  vhh = -17.248046773446934     (mV) 
  Ah = 1617.700731891178     (/ms) 
  b1h = -0.008950495371236813     (/mV) 
  c1h = -0.0002160257942825566     (/mV2) 
  d1h = 1.5760144473261236e-06     (/mV3) 
  b2h = 0.03184293962496223     (/mV) 
  c2h = -0.00017713711768543056     (/mV2) 
  d2h = -3.237252084317063e-06     (/mV3) 

  am = 0.06849296636369658     (/mV) 
  bm = -0.9794437916863961     (1) 
  vhm = -107.65472513187586     (mV) 
  Am = 35.14500237912589     (/ms) 
  b1m = -0.048842911677764854     (/mV) 
  c1m = 0.0004177714496913431     (/mV2) 
  d1m = -1.2015269858466243e-06     (/mV3) 
  b2m = -0.05682310713782679     (/mV) 
  c2m = 0.001345557972777253     (/mV2) 
  d2m = -5.8115129244024855e-06     (/mV3) 
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
  g = gbar*h*m*m
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