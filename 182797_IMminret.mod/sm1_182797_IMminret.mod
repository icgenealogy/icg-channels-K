NEURON
{
  SUFFIX IMminret 
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

  am = 0.142878665782301     (/mV) 
  bm = -8.573068475305128     (1) 
  vhm = -52.65481264647122     (mV) 
  Am = 270.71608338821005     (/ms) 
  b1m = 0.021898850706826182     (/mV) 
  c1m = -0.00047082025570705497     (/mV2) 
  d1m = 2.4198908246810383e-06     (/mV3) 
  b2m = 0.10468473742138765     (/mV) 
  c2m = -0.0009183285542410228     (/mV2) 
  d2m = 2.526271629772579e-06     (/mV3) 
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