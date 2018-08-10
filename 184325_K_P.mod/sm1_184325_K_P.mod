NEURON
{
  SUFFIX K_P 
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

  ah = -0.09087502154647435     (/mV) 
  bh = 4.908013408532564     (1) 
  vhh = -17.26414478244026     (mV) 
  Ah = 426.8626837928252     (/ms) 
  b1h = -0.008944405933943533     (/mV) 
  c1h = -0.00021612156405110167     (/mV2) 
  d1h = 1.5761472363486124e-06     (/mV3) 
  b2h = 0.03184111464925334     (/mV) 
  c2h = -0.00017679449717333833     (/mV2) 
  d2h = -3.2382921402260136e-06     (/mV3) 

  am = 0.06849296636369621     (/mV) 
  bm = -0.9794437916863571     (1) 
  vhm = -60.87747909862722     (mV) 
  Am = 20.20467281951215     (/ms) 
  b1m = 0.0948350581817801     (/mV) 
  c1m = 0.0028529330951639117     (/mV2) 
  d1m = 3.9319328879711544e-05     (/mV3) 
  b2m = 0.03376048764024993     (/mV) 
  c2m = -0.00011399135590803816     (/mV2) 
  d2m = 2.842775407319694e-07     (/mV3) 
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
  g = gbar*h*m*m
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