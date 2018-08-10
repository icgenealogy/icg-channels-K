NEURON
{
  SUFFIX imZ 
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

  am = 0.11112502464178127     (/mV) 
  bm = -3.334118210246568     (1) 
  vhm = -34.82833616161835     (mV) 
  Am = 338.3749690176624     (/ms) 
  b1m = -0.045194231952256296     (/mV) 
  c1m = 0.0003057818393148745     (/mV2) 
  d1m = -8.684150849277065e-07     (/mV3) 
  b2m = -0.06380060458125021     (/mV) 
  c2m = -0.0007854771085389497     (/mV2) 
  d2m = -4.392291024762103e-06     (/mV3) 
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