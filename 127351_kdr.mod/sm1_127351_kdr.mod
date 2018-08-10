NEURON
{
  SUFFIX kdr 
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

  an = 0.11222813215421547     (/mV) 
  bn = 1.4590350226381203     (1) 
  vhn = 15.882999853819591     (mV) 
  An = 46.5123149705672     (/ms) 
  b1n = -0.08595771637820575     (/mV) 
  c1n = -9.232470360376494e-05     (/mV2) 
  d1n = 7.035323320720395e-06     (/mV3) 
  b2n = -0.023329501884991773     (/mV) 
  c2n = 0.00019197961987783813     (/mV2) 
  d2n = 1.1082025852042745e-06     (/mV3) 
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