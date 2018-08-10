NEURON
{
  SUFFIX kx 
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

  am = 0.1881891186909873     (/mV) 
  bm = -8.112280803665955     (1) 
  vhm = -45.760174827095376     (mV) 
  Am = 105.98574047915537     (/ms) 
  b1m = -0.05086493189055587     (/mV) 
  c1m = 0.0003842919647139708     (/mV2) 
  d1m = -1.1635103613554398e-06     (/mV3) 
  b2m = -0.14023684915489776     (/mV) 
  c2m = -0.0023997830143114214     (/mV2) 
  d2m = -1.7394525779640173e-05     (/mV3) 
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