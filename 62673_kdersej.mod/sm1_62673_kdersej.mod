NEURON
{
  SUFFIX kder_sej 
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

  am = 0.11109795027788694     (/mV) 
  bm = -4.635261424602399     (1) 
  vhm = -40.33612822576944     (mV) 
  Am = 4.990202152905738     (/ms) 
  b1m = 0.03344265307129865     (/mV) 
  c1m = 0.0003331915489327361     (/mV2) 
  d1m = 1.8600546903778554e-06     (/mV3) 
  b2m = 0.07920728768611222     (/mV) 
  c2m = -0.0006504564155986731     (/mV2) 
  d2m = 2.0479790066928374e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
}

STATE
{
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*m*m*m
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}