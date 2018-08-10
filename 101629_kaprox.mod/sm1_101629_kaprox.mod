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

  al = -0.11223017693943452     (/mV) 
  bl = 3.5913596404421697     (1) 
  vhl = -11.971377244128705     (mV) 
  Al = 7.315386614700732     (/ms) 
  b1l = 0.01719374995085891     (/mV) 
  c1l = -5.1558207079258295e-06     (/mV2) 
  d1l = -2.6305317607930515e-07     (/mV3) 
  b2l = -0.07516430243092787     (/mV) 
  c2l = -0.0016961414798163021     (/mV2) 
  d2l = -5.870657957095992e-06     (/mV3) 

  an = 0.057987142491554095     (/mV) 
  bn = 2.0439159392680053     (1) 
  vhn = -17.32568048693127     (mV) 
  An = 0.816042821241069     (/ms) 
  b1n = 0.017721141322652796     (/mV) 
  c1n = -0.00016530849466739198     (/mV2) 
  d1n = -5.253075296063703e-07     (/mV3) 
  b2n = -0.0996373046925849     (/mV) 
  c2n = -0.0008381760126263447     (/mV2) 
  d2n = 1.107297421765285e-05     (/mV3) 
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