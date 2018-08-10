NEURON
{
  SUFFIX kacurrent 
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

  al = -0.1099990881421012     (/mV) 
  bl = 6.159939596493002     (1) 
  vhl = -35.92112876558254     (mV) 
  Al = 7.416245205560595     (/ms) 
  b1l = 0.07600259011307714     (/mV) 
  c1l = 0.0013100157485364007     (/mV2) 
  d1l = -1.0295401214981123e-05     (/mV3) 
  b2l = -0.01900783998622462     (/mV) 
  c2l = 5.4305735352585535e-05     (/mV2) 
  d2l = -3.017916391902783e-08     (/mV3) 

  an = 0.05831859707865206     (/mV) 
  bn = 0.6503931043866372     (1) 
  vhn = 29.979702649498947     (mV) 
  An = 3.3074183721081876     (/ms) 
  b1n = 0.008404030391495341     (/mV) 
  c1n = -0.00013875318452434608     (/mV2) 
  d1n = -3.59770349345432e-07     (/mV3) 
  b2n = 0.04226487099412021     (/mV) 
  c2n = -0.0001620170111511952     (/mV2) 
  d2n = 9.729500595416905e-07     (/mV3) 

  and = 0.07055708064363278     (/mV) 
  bnd = -0.0575231182109069     (1) 
  vhnd = -6.842396483868693     (mV) 
  And = 1.8996755668148713     (/ms) 
  b1nd = 0.04335117812256061     (/mV) 
  c1nd = -0.0001474688286473801     (/mV2) 
  d1nd = -3.0603333549041478e-06     (/mV3) 
  b2nd = 0.015857230600494748     (/mV) 
  c2nd = 0.00015670570486355444     (/mV2) 
  d2nd = -8.366313172808861e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  lInf 
  lTau 
  nInf 
  nTau 
  ndInf 
  ndTau 
}

STATE
{
  l
  n
  nd
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*l*l*n*n*nd
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  l' = (lInf - l) / lTau 
  n' = (nInf - n) / nTau 
  nd' = (ndInf - nd) / ndTau 
}

INITIAL
{
  rates(v)
  l = lInf 
  n = nInf 
  nd = ndInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    lInf = 1/(1 + exp(-al*v + bl)) 
    lTau = Al / ( exp(-(b1l*(v-vhl) + c1l*(v-vhl)^2 + d1l*(v-vhl)^3)) + exp((b2l*(v-vhl) + c2l*(v-vhl)^2 + d2l*(v-vhl)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    ndInf = 1/(1 + exp(-and*v + bnd)) 
    ndTau = And / ( exp(-(b1nd*(v-vhnd) + c1nd*(v-vhnd)^2 + d1nd*(v-vhnd)^3)) + exp((b2nd*(v-vhnd) + c2nd*(v-vhnd)^2 + d2nd*(v-vhnd)^3)) ) 


  UNITSON
}