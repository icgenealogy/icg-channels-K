NEURON
{
  SUFFIX KDRf 
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

  apfast = 0.11627595956946216     (/mV) 
  bpfast = -1.8836112131245593     (1) 
  vhpfast = -7.053139034803205     (mV) 
  Apfast = 7.121411015136257     (/ms) 
  b1pfast = -0.031196746618575435     (/mV) 
  c1pfast = -0.0004559115676798224     (/mV2) 
  d1pfast = 3.8007298166211013e-06     (/mV3) 
  b2pfast = -0.03949496400675493     (/mV) 
  c2pfast = 0.00036014769731448687     (/mV2) 
  d2pfast = 3.99355544081249e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  pfastInf 
  pfastTau 
}

STATE
{
  pfast
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*pfast
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  pfast' = (pfastInf - pfast) / pfastTau 
}

INITIAL
{
  rates(v)
  pfast = pfastInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    pfastInf = 1/(1 + exp(-apfast*v + bpfast)) 
    pfastTau = Apfast / ( exp(-(b1pfast*(v-vhpfast) + c1pfast*(v-vhpfast)^2 + d1pfast*(v-vhpfast)^3)) + exp((b2pfast*(v-vhpfast) + c2pfast*(v-vhpfast)^2 + d2pfast*(v-vhpfast)^3)) ) 


  UNITSON
}