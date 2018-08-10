NEURON
{
  SUFFIX Ksoma 
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

  an = 0.045145474264142386     (/mV) 
  bn = -0.7291819185089863     (1) 
  vhn = 16.56604915176502     (mV) 
  An = 1.0440208497673333     (/ms) 
  b1n = -0.0011202486874864596     (/mV) 
  c1n = -0.00011600702209187305     (/mV2) 
  d1n = -5.685772552428751e-07     (/mV3) 
  b2n = 0.02164931033827924     (/mV) 
  c2n = -1.2644333763112342e-05     (/mV2) 
  d2n = -1.1956297126023476e-06     (/mV3) 
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
  g = gbar*n*n*n*n
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