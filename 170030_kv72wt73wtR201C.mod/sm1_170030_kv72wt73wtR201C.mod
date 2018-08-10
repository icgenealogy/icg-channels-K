NEURON
{
  SUFFIX kv72wt73wtR201C 
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

  am = 0.07518757379850025     (/mV) 
  bm = -3.2631468587432186     (1) 
  vhm = -11.67512458531017     (mV) 
  Am = 15.802104941036658     (/ms) 
  b1m = 0.000949926391990859     (/mV) 
  c1m = -8.790505458579586e-05     (/mV2) 
  d1m = 6.814642359739548e-07     (/mV3) 
  b2m = -0.04965298441842528     (/mV) 
  c2m = -0.00428202989287739     (/mV2) 
  d2m = -4.212842955101975e-05     (/mV3) 
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
  g = gbar*m
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