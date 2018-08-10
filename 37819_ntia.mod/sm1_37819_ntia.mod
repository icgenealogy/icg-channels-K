NEURON
{
  SUFFIX iao 
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

  am = 0.13513486179595977     (/mV) 
  bm = -6.621604427621096     (1) 
  vhm = -76.76539174498772     (mV) 
  Am = 4.662265793467455     (/ms) 
  b1m = -0.05778850459905852     (/mV) 
  c1m = 0.00043860315274943433     (/mV2) 
  d1m = -1.0923203526780067e-06     (/mV3) 
  b2m = -0.07756994250506692     (/mV) 
  c2m = -0.0006419727658383782     (/mV2) 
  d2m = -3.3723434816887555e-06     (/mV3) 

  ah = -0.1999887420891276     (/mV) 
  bh = 18.599102601300736     (1) 
  vhh = -100.53846329587165     (mV) 
  Ah = 70.0069953467666     (/ms) 
  b1h = -3.103181501997206e-07     (/mV) 
  c1h = 1.1061793526790504e-08     (/mV2) 
  d1h = -4.411153714346147e-11     (/mV3) 
  b2h = -2.881327651485101     (/mV) 
  c2h = 0.005063659212720764     (/mV2) 
  d2h = -4.026671802705155e-05     (/mV3) 

  an = -0.1999875526655743     (/mV) 
  bn = 18.598927893412146     (1) 
  vhn = -101.9816885564249     (mV) 
  An = 62.885028298139375     (/ms) 
  b1n = -0.0014366527694745638     (/mV) 
  c1n = 1.2761112594014136e-05     (/mV2) 
  d1n = -3.428716443576808e-08     (/mV3) 
  b2n = -1.4355807387144808     (/mV) 
  c2n = 0.10315557554071557     (/mV2) 
  d2n = -0.002062351900248949     (/mV3) 
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
  hInf 
  hTau 
  nInf 
  nTau 
}

STATE
{
  m
  h
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*m*m*m*m*m*h*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  m = mInf 
  h = hInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}