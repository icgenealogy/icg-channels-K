NEURON
{
  SUFFIX kap 
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

  al = -0.11223020247549506     (/mV) 
  bl = 5.611507436747995     (1) 
  vhl = -29.477837050762428     (mV) 
  Al = 7.619625545070251     (/ms) 
  b1l = 0.07621052372898791     (/mV) 
  c1l = 0.0011929128335295344     (/mV2) 
  d1l = -9.675497803122447e-06     (/mV3) 
  b2l = -0.01821783662048014     (/mV) 
  c2l = 4.425818434440342e-05     (/mV2) 
  d2l = 1.3289827921016313e-08     (/mV3) 

  an = 0.05796171833635216     (/mV) 
  bn = 0.9998513392145176     (1) 
  vhn = -35.57215850037624     (mV) 
  An = 0.7903185282228372     (/ms) 
  b1n = 0.09852852645719737     (/mV) 
  c1n = 0.0010224965696954782     (/mV2) 
  d1n = -1.1591826307564747e-05     (/mV3) 
  b2n = -0.01993643916776327     (/mV) 
  c2n = 0.00021836666451694096     (/mV2) 
  d2n = 1.4435407439760777e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  lInf 
  lTau 
  nInf 
  nTau 
}

STATE
{
  l
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*l*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  l' = (lInf - l) / lTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  l = lInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    lInf = 1/(1 + exp(-al*v + bl)) 
    lTau = Al / ( exp(-(b1l*(v-vhl) + c1l*(v-vhl)^2 + d1l*(v-vhl)^3)) + exp((b2l*(v-vhl) + c2l*(v-vhl)^2 + d2l*(v-vhl)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}