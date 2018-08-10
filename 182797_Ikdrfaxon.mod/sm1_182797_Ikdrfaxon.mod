NEURON
{
  SUFFIX Ikdrfaxon 
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

  ah = -0.10147297097797246     (/mV) 
  bh = 3.8874819904987943     (1) 
  vhh = -43.82093206351344     (mV) 
  Ah = 2000.0016282415017     (/ms) 
  b1h = -0.00015565763541947766     (/mV) 
  c1h = 4.974006870686776e-07     (/mV2) 
  d1h = -2.8930259218152974e-11     (/mV3) 
  b2h = -0.00015538094275312025     (/mV) 
  c2h = 4.6921643689608045e-07     (/mV2) 
  d2h = 9.736318640140907e-11     (/mV3) 

  am = 0.08149469725653388     (/mV) 
  bm = -0.6667898124617956     (1) 
  vhm = -33.765563180919415     (mV) 
  Am = 5.884186582229372     (/ms) 
  b1m = -0.028371492963255234     (/mV) 
  c1m = -2.2603046228961953e-05     (/mV2) 
  d1m = 1.1635854524014314e-07     (/mV3) 
  b2m = -0.07215936640121344     (/mV) 
  c2m = -4.9534542990266266e-05     (/mV2) 
  d2m = -1.0450229353252068e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hInf 
  hTau 
  mInf 
  mTau 
}

STATE
{
  h
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*m
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}