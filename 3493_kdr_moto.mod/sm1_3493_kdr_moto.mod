NEURON
{
  SUFFIX kdrmot 
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

  am = 0.1176425027539367     (/mV) 
  bm = -5.152677361215342     (1) 
  vhm = -50.19253710156249     (mV) 
  Am = 24.184676006969994     (/ms) 
  b1m = 0.05693834267492117     (/mV) 
  c1m = 0.0001408752748590369     (/mV2) 
  d1m = 1.138420384861555e-06     (/mV3) 
  b2m = 0.054487819514533825     (/mV) 
  c2m = 0.0001427614548649132     (/mV2) 
  d2m = -9.90757588780004e-07     (/mV3) 
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