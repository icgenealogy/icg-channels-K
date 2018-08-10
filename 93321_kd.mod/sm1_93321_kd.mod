NEURON
{
  SUFFIX kd 
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

  am = 0.08460564699693683     (/mV) 
  bm = -1.0391271812427754     (1) 
  vhm = -89.01833778935139     (mV) 
  Am = 14.122544230844563     (/ms) 
  b1m = 0.0022421584616576537     (/mV) 
  c1m = -0.00026545445914227876     (/mV2) 
  d1m = 1.0252878179582246e-06     (/mV3) 
  b2m = -0.0022395507154303602     (/mV) 
  c2m = 0.0002654742857299946     (/mV2) 
  d2m = -1.0253626764246554e-06     (/mV3) 
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