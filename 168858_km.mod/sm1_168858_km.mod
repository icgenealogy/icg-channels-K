NEURON
{
  SUFFIX km 
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

  am = 0.1936026322857213     (/mV) 
  bm = -5.078801796467555     (1) 
  vhm = -20.549703055032328     (mV) 
  Am = 179.6918659769188     (/ms) 
  b1m = -0.10410128118908629     (/mV) 
  c1m = 0.0027372146617201825     (/mV2) 
  d1m = -1.592218010743913e-05     (/mV3) 
  b2m = -0.007269595045353455     (/mV) 
  c2m = 0.0006331624933113009     (/mV2) 
  d2m = -4.388505738571426e-06     (/mV3) 
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