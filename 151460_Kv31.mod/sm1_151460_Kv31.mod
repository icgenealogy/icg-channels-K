NEURON
{
  SUFFIX Kv31 
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

  ap = 0.11110298792260061     (/mV) 
  bp = 0.009456397011601347     (1) 
  vhp = -15.109672863098302     (mV) 
  Ap = 16.87719109451895     (/ms) 
  b1p = -0.07846654166961825     (/mV) 
  c1p = 0.0006400633592534586     (/mV2) 
  d1p = -1.5702507292386862e-06     (/mV3) 
  b2p = -0.05694185849727425     (/mV) 
  c2p = 8.019468020520303e-05     (/mV2) 
  d2p = 4.0922435244715435e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  pInf 
  pTau 
}

STATE
{
  p
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*p
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  p' = (pInf - p) / pTau 
}

INITIAL
{
  rates(v)
  p = pInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    pInf = 1/(1 + exp(-ap*v + bp)) 
    pTau = Ap / ( exp(-(b1p*(v-vhp) + c1p*(v-vhp)^2 + d1p*(v-vhp)^3)) + exp((b2p*(v-vhp) + c2p*(v-vhp)^2 + d2p*(v-vhp)^3)) ) 


  UNITSON
}