NEURON
{
  SUFFIX km_hu 
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

  am = 0.13094153401917336     (/mV) 
  bm = -5.63058431282904     (1) 
  vhm = -42.227562504227144     (mV) 
  Am = 112.47822667371067     (/ms) 
  b1m = -0.06574243608730258     (/mV) 
  c1m = -6.893213249685954e-05     (/mV2) 
  d1m = 1.9740118475133527e-06     (/mV3) 
  b2m = -0.05935320442200959     (/mV) 
  c2m = 0.00021403732682331782     (/mV2) 
  d2m = 3.0007469230744538e-06     (/mV3) 
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