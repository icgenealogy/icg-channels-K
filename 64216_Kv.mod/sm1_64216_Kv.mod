NEURON
{
  SUFFIX Kv 
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

  ahKv = -0.1277713164608074     (/mV) 
  bhKv = -0.669187029332726     (1) 
  vhhKv = 24.060992125512474     (mV) 
  AhKv = 4798.976489071848     (/ms) 
  b1hKv = -0.028560117356898167     (/mV) 
  c1hKv = 0.002140669762857535     (/mV2) 
  d1hKv = -2.0612933611365048e-05     (/mV3) 
  b2hKv = 0.006010316829371231     (/mV) 
  c2hKv = 0.00025929980211251866     (/mV2) 
  d2hKv = -3.261492903398366e-06     (/mV3) 

  amKv = 0.04207355660012541     (/mV) 
  bmKv = -1.2323315226199878     (1) 
  vhmKv = -22.85886944775399     (mV) 
  AmKv = 32.73920382462793     (/ms) 
  b1mKv = 0.021437026518956362     (/mV) 
  c1mKv = -3.187610984029264e-05     (/mV2) 
  d1mKv = -1.0920998574726457e-07     (/mV3) 
  b2mKv = 0.019326518600527497     (/mV) 
  c2mKv = -3.3198095442697504e-05     (/mV2) 
  d2mKv = 2.057547888411247e-08     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hKvInf 
  hKvTau 
  mKvInf 
  mKvTau 
}

STATE
{
  hKv
  mKv
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*hKv*mKv*mKv*mKv
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  hKv' = (hKvInf - hKv) / hKvTau 
  mKv' = (mKvInf - mKv) / mKvTau 
}

INITIAL
{
  rates(v)
  hKv = hKvInf 
  mKv = mKvInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hKvInf = 1/(1 + exp(-ahKv*v + bhKv)) 
    hKvTau = AhKv / ( exp(-(b1hKv*(v-vhhKv) + c1hKv*(v-vhhKv)^2 + d1hKv*(v-vhhKv)^3)) + exp((b2hKv*(v-vhhKv) + c2hKv*(v-vhhKv)^2 + d2hKv*(v-vhhKv)^3)) ) 

    mKvInf = 1/(1 + exp(-amKv*v + bmKv)) 
    mKvTau = AmKv / ( exp(-(b1mKv*(v-vhmKv) + c1mKv*(v-vhmKv)^2 + d1mKv*(v-vhmKv)^3)) + exp((b2mKv*(v-vhmKv) + c2mKv*(v-vhmKv)^2 + d2mKv*(v-vhmKv)^3)) ) 


  UNITSON
}