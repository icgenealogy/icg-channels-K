NEURON
{
  SUFFIX Ito 
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

  an = -0.1886552604649436     (/mV) 
  bn = 8.131210730737417     (1) 
  vhn = -50.895440231591174     (mV) 
  An = 159.4812470239739     (/ms) 
  b1n = 0.14464291016434977     (/mV) 
  c1n = 0.002860236206233903     (/mV2) 
  d1n = 2.5535160636754948e-05     (/mV3) 
  b2n = 0.10360322741428699     (/mV) 
  c2n = -0.0012055056324820604     (/mV2) 
  d2n = 4.2775636725862974e-06     (/mV3) 

  am = 0.057012459219310634     (/mV) 
  bm = -1.1670422251223058     (1) 
  vhm = -36.38636911363178     (mV) 
  Am = 10.480144373438213     (/ms) 
  b1m = -0.09888083653935038     (/mV) 
  c1m = 0.0011054960568508889     (/mV2) 
  d1m = -4.60179705015195e-06     (/mV3) 
  b2m = -0.06964209609015615     (/mV) 
  c2m = -0.0008855817104022496     (/mV2) 
  d2m = -4.392664918489138e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  nInf 
  nTau 
  mInf 
  mTau 
}

STATE
{
  n
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n*m*m*m
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  n = nInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}