NEURON
{
  SUFFIX kdrRL 
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

  am = 0.049994709260362966     (/mV) 
  bm = -1.2497770336545795     (1) 
  vhm = -45.089643244108224     (mV) 
  Am = 7.439660361310075     (/ms) 
  b1m = 0.1619034905898897     (/mV) 
  c1m = 0.004893034723144463     (/mV2) 
  d1m = 4.6102237339796606e-05     (/mV3) 
  b2m = 0.06996398679579681     (/mV) 
  c2m = -0.0008515219765356912     (/mV2) 
  d2m = 3.150508119402917e-06     (/mV3) 
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