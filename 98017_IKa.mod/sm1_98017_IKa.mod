NEURON
{
  SUFFIX IKa 
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

  aa = 0.11764614990467921     (/mV) 
  ba = -7.058772899651177     (1) 
  vha = -66.99897062960326     (mV) 
  Aa = 79.07081427687434     (/ms) 
  b1a = -0.23183326141371322     (/mV) 
  c1a = 0.0025684279424741207     (/mV2) 
  d1a = -8.326889428823558e-06     (/mV3) 
  b2a = -0.012320144382881038     (/mV) 
  c2a = 0.00025392496316931567     (/mV2) 
  d2a = -1.1082688662773481e-07     (/mV3) 

  ab = -0.16665134865585773     (/mV) 
  bb = 12.998978373433117     (1) 
  vhb = -64.51924214385568     (mV) 
  Ab = 100.57246449067551     (/ms) 
  b1b = -0.10644604375524065     (/mV) 
  c1b = 0.0019861681393613735     (/mV2) 
  d1b = -8.426752547185665e-06     (/mV3) 
  b2b = -0.0010505564119421124     (/mV) 
  c2b = 0.0003702598282809043     (/mV2) 
  d2b = -2.00578194301889e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  aInf 
  aTau 
  bInf 
  bTau 
}

STATE
{
  a
  b
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*a*b
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  a' = (aInf - a) / aTau 
  b' = (bInf - b) / bTau 
}

INITIAL
{
  rates(v)
  a = aInf 
  b = bInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    aInf = 1/(1 + exp(-aa*v + ba)) 
    aTau = Aa / ( exp(-(b1a*(v-vha) + c1a*(v-vha)^2 + d1a*(v-vha)^3)) + exp((b2a*(v-vha) + c2a*(v-vha)^2 + d2a*(v-vha)^3)) ) 

    bInf = 1/(1 + exp(-ab*v + bb)) 
    bTau = Ab / ( exp(-(b1b*(v-vhb) + c1b*(v-vhb)^2 + d1b*(v-vhb)^3)) + exp((b2b*(v-vhb) + c2b*(v-vhb)^2 + d2b*(v-vhb)^3)) ) 


  UNITSON
}