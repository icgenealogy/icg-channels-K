NEURON
{
  SUFFIX KDR 
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

  an = 0.10574559678084722     (/mV) 
  bn = -4.043397680477977     (1) 
  vhn = -36.71169049701654     (mV) 
  An = 0.397397828816839     (/ms) 
  b1n = -0.07280127194203016     (/mV) 
  c1n = 0.0003748681470162273     (/mV2) 
  d1n = 1.0342432546880925e-06     (/mV3) 
  b2n = -0.022834315715882805     (/mV) 
  c2n = 0.0003286249796786395     (/mV2) 
  d2n = 3.5308423688096365e-06     (/mV3) 
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