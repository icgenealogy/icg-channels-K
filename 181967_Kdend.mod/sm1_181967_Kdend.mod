NEURON
{
  SUFFIX Kdend 
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

  an = 0.05184634242718803     (/mV) 
  bn = -0.6935143277972614     (1) 
  vhn = -1.5808741968062352     (mV) 
  An = 5.901579614554657     (/ms) 
  b1n = -0.02281731531950304     (/mV) 
  c1n = -6.19892238223776e-05     (/mV2) 
  d1n = 1.3369085954134277e-06     (/mV3) 
  b2n = -0.005748289911499877     (/mV) 
  c2n = 0.00010389145505215244     (/mV2) 
  d2n = 8.661947295803605e-07     (/mV3) 
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