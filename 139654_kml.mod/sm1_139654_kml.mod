NEURON
{
  SUFFIX kml 
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

  an = 0.19915244782226235     (/mV) 
  bn = -3.9737283459076584     (1) 
  vhn = -18.792817243059226     (mV) 
  An = 5.981836014299206     (/ms) 
  b1n = 0.04607766359382587     (/mV) 
  c1n = -8.23852383353084e-05     (/mV2) 
  d1n = -6.443457522478314e-07     (/mV3) 
  b2n = 0.05148587911938405     (/mV) 
  c2n = -1.0214783325811887e-05     (/mV2) 
  d2n = -1.4117677000557033e-07     (/mV3) 
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