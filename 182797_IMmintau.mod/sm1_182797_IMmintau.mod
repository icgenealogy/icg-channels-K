NEURON
{
  SUFFIX IM 
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

  am = 0.1428547796512211     (/mV) 
  bm = -3.857107818311156     (1) 
  vhm = -52.62451658396347     (mV) 
  Am = 270.0734273338434     (/ms) 
  b1m = -0.10087543392035701     (/mV) 
  c1m = 0.0007432034918958934     (/mV2) 
  d1m = -9.730796702231523e-07     (/mV3) 
  b2m = -0.01973672639052958     (/mV) 
  c2m = 0.0005543870247605471     (/mV2) 
  d2m = -1.7419621153440547e-06     (/mV3) 
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