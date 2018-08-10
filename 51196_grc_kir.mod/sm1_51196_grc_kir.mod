NEURON
{
  SUFFIX GrC_Kir 
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

  ad = -0.06900012225811676     (/mV) 
  bd = 6.037796426964715     (1) 
  vhd = -96.72120569882442     (mV) 
  Ad = 0.9139723229552105     (/ms) 
  b1d = 0.05043846133827987     (/mV) 
  c1d = 0.0003030915554859871     (/mV2) 
  d1d = -3.240124480788317e-06     (/mV3) 
  b2d = 0.01624000302588027     (/mV) 
  c2d = 0.00012303024042864714     (/mV2) 
  d2d = -5.279504645212388e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  dInf 
  dTau 
}

STATE
{
  d
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*d
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  d' = (dInf - d) / dTau 
}

INITIAL
{
  rates(v)
  d = dInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    dInf = 1/(1 + exp(-ad*v + bd)) 
    dTau = Ad / ( exp(-(b1d*(v-vhd) + c1d*(v-vhd)^2 + d1d*(v-vhd)^3)) + exp((b2d*(v-vhd) + c2d*(v-vhd)^2 + d2d*(v-vhd)^3)) ) 


  UNITSON
}