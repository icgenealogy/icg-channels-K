NEURON
{
  SUFFIX ka_chan 
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

  ah = -0.19999896197880584     (/mV) 
  bh = 15.999926691948092     (1) 
  vhh = -51.69508268630331     (mV) 
  Ah = 12.02000001048143     (/ms) 
  b1h = 1.8476056690707477e-09     (/mV) 
  c1h = -6.681048077042069e-09     (/mV2) 
  d1h = -3.437245525129847e-11     (/mV3) 
  b2h = 1.909124839782316e-09     (/mV) 
  c2h = -6.6791932610206276e-09     (/mV2) 
  d2h = -3.4406747820854235e-11     (/mV3) 

  am = 0.1299996805380268     (/mV) 
  bm = -9.099978292569673     (1) 
  vhm = -43.312794023895066     (mV) 
  Am = 1.9728220128459135     (/ms) 
  b1m = -0.06294907135801567     (/mV) 
  c1m = 0.0005512514002605782     (/mV2) 
  d1m = -1.6250142623362555e-06     (/mV3) 
  b2m = 0.012075524626081607     (/mV) 
  c2m = 0.0005842739246612964     (/mV2) 
  d2m = -7.745145186205477e-06     (/mV3) 
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
  g = gbar*h*m*m*m
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