NEURON
{
  SUFFIX Ks 
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

  aa = 0.15384598376212455     (/mV) 
  ba = -5.230755780226697     (1) 
  vha = -51.69508268630331     (mV) 
  Aa = 12.02000001048143     (/ms) 
  b1a = 1.8476056690707477e-09     (/mV) 
  c1a = -6.681048077042069e-09     (/mV2) 
  d1a = -3.437245525129847e-11     (/mV3) 
  b2a = 1.909124839782316e-09     (/mV) 
  c2a = -6.6791932610206276e-09     (/mV2) 
  d2a = -3.4406747820854235e-11     (/mV3) 

  ab = -0.15144041915392334     (/mV) 
  bb = 9.844533425761618     (1) 
  vhb = -26.4373340947171     (mV) 
  Ab = 852.6325650266455     (/ms) 
  b1b = 0.000565093392043079     (/mV) 
  c1b = 3.4570482911653265e-06     (/mV2) 
  d1b = 1.7023369224614825e-06     (/mV3) 
  b2b = -0.0005192300812407747     (/mV) 
  c2b = 0.00014353783132176268     (/mV2) 
  d2b = -7.618609152396543e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  aInf 
  aTau 
  bInf 
  bTau 
}

STATE
{
  a
  b
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*a*b
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  a' = (aInf - a) / aTau 
  b' = (bInf - b) / bTau 
}

INITIAL
{
  rates(v)
  a = aInf 
  b = bInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    aInf = 1/(1 + exp(-aa*v + ba)) 
    aTau = Aa / ( exp(-(b1a*(v-vha) + c1a*(v-vha)^2 + d1a*(v-vha)^3)) + exp((b2a*(v-vha) + c2a*(v-vha)^2 + d2a*(v-vha)^3)) ) 

    bInf = 1/(1 + exp(-ab*v + bb)) 
    bTau = Ab / ( exp(-(b1b*(v-vhb) + c1b*(v-vhb)^2 + d1b*(v-vhb)^3)) + exp((b2b*(v-vhb) + c2b*(v-vhb)^2 + d2b*(v-vhb)^3)) ) 


  UNITSON
}