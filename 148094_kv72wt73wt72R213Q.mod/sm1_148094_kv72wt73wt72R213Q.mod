NEURON
{
  SUFFIX kvR213Q 
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

  am = 0.06872647465568521     (/mV) 
  bm = -1.009568444554577     (1) 
  vhm = -8.353272648685769     (mV) 
  Am = 24.63207471466056     (/ms) 
  b1m = 0.10398382494265625     (/mV) 
  c1m = 0.0009583270143214597     (/mV2) 
  d1m = -1.692237427420952e-05     (/mV3) 
  b2m = 0.0010647519241293486     (/mV) 
  c2m = 0.0002621268252535673     (/mV2) 
  d2m = -2.2151624345433537e-06     (/mV3) 
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