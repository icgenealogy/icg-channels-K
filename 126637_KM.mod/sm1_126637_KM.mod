NEURON
{
  SUFFIX KM 
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

  am = 0.10000281672448384     (/mV) 
  bm = -3.5001495111510748     (1) 
  vhm = -34.90799630761442     (mV) 
  Am = 60.685159727418814     (/ms) 
  b1m = -0.025126684585324546     (/mV) 
  c1m = 1.6566319623209442e-06     (/mV2) 
  d1m = -6.230681901846636e-09     (/mV3) 
  b2m = -0.049863978673724     (/mV) 
  c2m = 1.90782730101103e-06     (/mV2) 
  d2m = 1.9387322455356264e-08     (/mV3) 
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