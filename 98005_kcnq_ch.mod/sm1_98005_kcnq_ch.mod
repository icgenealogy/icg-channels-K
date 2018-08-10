NEURON
{
  SUFFIX kcnq_ch 
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

  ac = -0.10286395549171559     (/mV) 
  bc = 5.045068799314814     (1) 
  vhc = -11.521568580608767     (mV) 
  Ac = 21.71051807051059     (/ms) 
  b1c = -0.0037965251521677214     (/mV) 
  c1c = -3.4532741014943546e-05     (/mV2) 
  d1c = 2.781886127397647e-07     (/mV3) 
  b2c = 0.003797933268951972     (/mV) 
  c2c = -0.0002689385147235329     (/mV2) 
  d2c = 1.559472673828969e-06     (/mV3) 

  ao = 0.10286395548704634     (/mV) 
  bo = -5.04506879905542     (1) 
  vho = -11.521568580608767     (mV) 
  Ao = 21.71051807051059     (/ms) 
  b1o = -0.0037965251521677214     (/mV) 
  c1o = -3.4532741014943546e-05     (/mV2) 
  d1o = 2.781886127397647e-07     (/mV3) 
  b2o = 0.003797933268951972     (/mV) 
  c2o = -0.0002689385147235329     (/mV2) 
  d2o = 1.559472673828969e-06     (/mV3) 
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