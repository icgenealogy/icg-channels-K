NEURON
{
  SUFFIX kdrsoma 
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

  am = 0.0998287056198459     (/mV) 
  bm = -1.9450775813771495     (1) 
  vhm = -10.411073527058361     (mV) 
  Am = 8.335017452987032     (/ms) 
  b1m = -0.16417269087109035     (/mV) 
  c1m = 0.002496516824528124     (/mV2) 
  d1m = -1.2022871269045605e-05     (/mV3) 
  b2m = -0.182778425872854     (/mV) 
  c2m = -0.0033059808044028317     (/mV2) 
  d2m = -1.9539585324272616e-05     (/mV3) 
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