NEURON
{
  SUFFIX kdrmt 
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

  am = 0.09999906955319282     (/mV) 
  bm = 2.100016748226314     (1) 
  vhm = -49.240481262698985     (mV) 
  Am = 68.48667066694584     (/ms) 
  b1m = 0.02687118134957374     (/mV) 
  c1m = -6.652615534004444e-06     (/mV2) 
  d1m = -1.952450174096538e-08     (/mV3) 
  b2m = 0.02801508695646497     (/mV) 
  c2m = -4.602063844144461e-06     (/mV2) 
  d2m = 1.3301674850143907e-08     (/mV3) 
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