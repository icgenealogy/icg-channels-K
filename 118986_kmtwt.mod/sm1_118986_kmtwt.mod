NEURON
{
  SUFFIX kmtwt 
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

  am = 0.08196665154119497     (/mV) 
  bm = -2.655722427141616     (1) 
  vhm = -56.221433488938956     (mV) 
  Am = 40.61392252257766     (/ms) 
  b1m = 0.12627381719481134     (/mV) 
  c1m = 0.00322387243831608     (/mV2) 
  d1m = 3.513556621508876e-05     (/mV3) 
  b2m = 0.07795108395342763     (/mV) 
  c2m = -0.0008985131717757473     (/mV2) 
  d2m = 3.1558211506800114e-06     (/mV3) 
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