NEURON
{
  SUFFIX IKdrSM 
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

  am = 0.08256924478375016     (/mV) 
  bm = -2.628831725249166     (1) 
  vhm = -43.0732759962157     (mV) 
  Am = 15.703507323547937     (/ms) 
  b1m = -0.041957261704528244     (/mV) 
  c1m = 0.0002504642224201161     (/mV2) 
  d1m = -6.345848227150597e-07     (/mV3) 
  b2m = -0.05308169105501327     (/mV) 
  c2m = -0.0004428236811007563     (/mV2) 
  d2m = -3.287058880384352e-06     (/mV3) 
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