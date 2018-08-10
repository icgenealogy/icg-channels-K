NEURON
{
  SUFFIX kavn 
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

  ah = -0.01666467193867362     (/mV) 
  bh = 1.9394794831051824e-05     (1) 
  vhh = -13.175717370828204     (mV) 
  Ah = 158.25173477039183     (/ms) 
  b1h = 0.00023141169085946534     (/mV) 
  c1h = -0.0001711842639588113     (/mV2) 
  d1h = 8.949999753783503e-07     (/mV3) 
  b2h = -0.00023418865428403323     (/mV) 
  c2h = -1.9438980445762988e-05     (/mV2) 
  d2h = -1.252734124412483e-06     (/mV3) 

  am = 0.07812468634326214     (/mV) 
  bm = -1.6406090866814704     (1) 
  vhm = -75.01858787373455     (mV) 
  Am = 54.932542425707815     (/ms) 
  b1m = 0.1642203809670889     (/mV) 
  c1m = 1.1113990175052782e-05     (/mV2) 
  d1m = -1.35089992057459e-07     (/mV3) 
  b2m = 0.03587511590326662     (/mV) 
  c2m = 2.480338130499998e-06     (/mV2) 
  d2m = -1.7616335926017242e-08     (/mV3) 
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