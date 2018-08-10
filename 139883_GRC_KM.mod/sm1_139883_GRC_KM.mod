NEURON
{
  SUFFIX GRC_KM 
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

  an = 0.16667436338324798     (/mV) 
  bn = -5.833649873320012     (1) 
  vhn = 0.7425735761133218     (mV) 
  An = 48.46511558934131     (/ms) 
  b1n = -0.04752713389796592     (/mV) 
  c1n = 0.0006297127805221387     (/mV2) 
  d1n = -1.5738315870847217e-06     (/mV3) 
  b2n = -0.006639789892931372     (/mV) 
  c2n = 0.00023042931377651195     (/mV2) 
  d2n = -5.880228285255849e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  nInf 
  nTau 
}

STATE
{
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}