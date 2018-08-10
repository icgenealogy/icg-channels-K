NEURON
{
  SUFFIX fKdr 
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

  am = 0.12820490153401545     (/mV) 
  bm = -5.128194835441416     (1) 
  vhm = -38.48930113707086     (mV) 
  Am = 13.97408773464591     (/ms) 
  b1m = -0.09256680438490664     (/mV) 
  c1m = 0.0002516554756230524     (/mV2) 
  d1m = 1.5845920239951816e-06     (/mV3) 
  b2m = -0.06747044592632796     (/mV) 
  c2m = 0.00022717715007479437     (/mV2) 
  d2m = 3.7223501970344884e-06     (/mV3) 
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