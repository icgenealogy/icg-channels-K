NEURON
{
  SUFFIX kdr_fs 
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

  am = 0.08682054688629917     (/mV) 
  bm = -2.342863091865878     (1) 
  vhm = -10.375192430215291     (mV) 
  Am = 8.292522104086068     (/ms) 
  b1m = -0.16375709279328418     (/mV) 
  c1m = 0.0024957059043605296     (/mV2) 
  d1m = -1.1961830511860753e-05     (/mV3) 
  b2m = -0.18148419486246528     (/mV) 
  c2m = -0.00326554952751493     (/mV2) 
  d2m = -1.9220082579062517e-05     (/mV3) 
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
  g = gbar*m*m*m*m
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