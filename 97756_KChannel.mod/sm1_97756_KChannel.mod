NEURON
{
  SUFFIX KCHANNEL 
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

  ac = -0.029999383312054664     (/mV) 
  bc = -0.6931698616175053     (1) 
  vhc = 82.0871628689181     (mV) 
  Ac = 75.21138638613823     (/ms) 
  b1c = -0.008062320782068885     (/mV) 
  c1c = 4.221282289454594e-05     (/mV2) 
  d1c = 5.812169799396396e-07     (/mV3) 
  b2c = 0.0031632623607146136     (/mV) 
  c2c = 0.0001468197453866051     (/mV2) 
  d2c = 3.382156347249516e-07     (/mV3) 

  ao = 0.02999938331204684     (/mV) 
  bo = 0.6931698616170198     (1) 
  vho = 82.0871628689181     (mV) 
  Ao = 75.21138638613823     (/ms) 
  b1o = -0.008062320782068885     (/mV) 
  c1o = 4.221282289454594e-05     (/mV2) 
  d1o = 5.812169799396396e-07     (/mV3) 
  b2o = 0.0031632623607146136     (/mV) 
  c2o = 0.0001468197453866051     (/mV2) 
  d2o = 3.382156347249516e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  cInf 
  cTau 
  oInf 
  oTau 
}

STATE
{
  c
  o
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*o
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  c' = (cInf - c) / cTau 
  o' = (oInf - o) / oTau 
}

INITIAL
{
  rates(v)
  c = cInf 
  o = oInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    cInf = 1/(1 + exp(-ac*v + bc)) 
    cTau = Ac / ( exp(-(b1c*(v-vhc) + c1c*(v-vhc)^2 + d1c*(v-vhc)^3)) + exp((b2c*(v-vhc) + c2c*(v-vhc)^2 + d2c*(v-vhc)^3)) ) 

    oInf = 1/(1 + exp(-ao*v + bo)) 
    oTau = Ao / ( exp(-(b1o*(v-vho) + c1o*(v-vho)^2 + d1o*(v-vho)^3)) + exp((b2o*(v-vho) + c2o*(v-vho)^2 + d2o*(v-vho)^3)) ) 


  UNITSON
}