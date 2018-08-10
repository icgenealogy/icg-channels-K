NEURON
{
  SUFFIX Kv 
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

  amKv = 0.03810049899040846     (/mV) 
  bmKv = -0.007917062756772793     (1) 
  vhmKv = -40.043138370109396     (mV) 
  AmKv = 16.145448584994753     (/ms) 
  b1mKv = 0.015272068567409697     (/mV) 
  c1mKv = 2.054602233627568e-06     (/mV2) 
  d1mKv = -1.0019641877193077e-06     (/mV3) 
  b2mKv = 0.00529640583785833     (/mV) 
  c2mKv = 0.00016274039818056348     (/mV2) 
  d2mKv = -8.341416000786875e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mKvInf 
  mKvTau 
}

STATE
{
  mKv
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*mKv*mKv*mKv*mKv
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  mKv' = (mKvInf - mKv) / mKvTau 
}

INITIAL
{
  rates(v)
  mKv = mKvInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mKvInf = 1/(1 + exp(-amKv*v + bmKv)) 
    mKvTau = AmKv / ( exp(-(b1mKv*(v-vhmKv) + c1mKv*(v-vhmKv)^2 + d1mKv*(v-vhmKv)^3)) + exp((b2mKv*(v-vhmKv) + c2mKv*(v-vhmKv)^2 + d2mKv*(v-vhmKv)^3)) ) 


  UNITSON
}