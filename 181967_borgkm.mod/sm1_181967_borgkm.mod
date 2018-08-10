NEURON
{
  SUFFIX borgkm 
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

  am = 0.37410297462438696     (/mV) 
  bm = -20.57564767490406     (1) 
  vhm = -55.00492010616366     (mV) 
  Am = 17.50623355156917     (/ms) 
  b1m = -0.022412489159952042     (/mV) 
  c1m = -2.2513779519564493e-07     (/mV2) 
  d1m = 1.726933084349051e-09     (/mV3) 
  b2m = -0.35193378437258055     (/mV) 
  c2m = -9.007787907283305e-05     (/mV2) 
  d2m = -1.917655822673749e-06     (/mV3) 
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