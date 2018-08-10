NEURON
{
  SUFFIX kmtquad 
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

  am = 0.07633519392389929     (/mV) 
  bm = -2.1144970411004453     (1) 
  vhm = -26.239087213410702     (mV) 
  Am = 32.26697048713615     (/ms) 
  b1m = -0.004481266302917672     (/mV) 
  c1m = -0.000311226467550466     (/mV2) 
  d1m = 2.2830473574271186e-06     (/mV3) 
  b2m = 0.061191816558708555     (/mV) 
  c2m = -0.0015475237397391563     (/mV2) 
  d2m = 8.625432835936676e-06     (/mV3) 
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