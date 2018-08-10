NEURON
{
  SUFFIX Im 
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

  am = 0.2000057280323602     (/mV) 
  bm = -7.000261092291547     (1) 
  vhm = -35.01246954747373     (mV) 
  Am = 79.94236652987672     (/ms) 
  b1m = -0.09973655687881677     (/mV) 
  c1m = -1.51577028782413e-05     (/mV2) 
  d1m = 2.8473100768039794e-07     (/mV3) 
  b2m = -0.09998314876681885     (/mV) 
  c2m = 6.394042452783638e-06     (/mV2) 
  d2m = 1.7901858097720392e-07     (/mV3) 
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