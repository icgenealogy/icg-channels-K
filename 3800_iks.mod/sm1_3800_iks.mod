NEURON
{
  SUFFIX IKs 
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

  am = 0.05795446764276315     (/mV) 
  bm = 0.2518897115323701     (1) 
  vhm = 45.56205595985307     (mV) 
  Am = 736.9822043922314     (/ms) 
  b1m = -0.0013155116688416298     (/mV) 
  c1m = -0.00030868926476894707     (/mV2) 
  d1m = -1.5773921522016799e-06     (/mV3) 
  b2m = 0.037341849056903895     (/mV) 
  c2m = -0.0003353170906057635     (/mV2) 
  d2m = -4.304449143197948e-06     (/mV3) 
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
  g = gbar*m*m
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