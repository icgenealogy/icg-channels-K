NEURON
{
  SUFFIX KDRI 
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

  ah = -0.01586343331253105     (/mV) 
  bh = -0.9297124730825106     (1) 
  vhh = -65.8856379138918     (mV) 
  Ah = 4.140000002255323     (/ms) 
  b1h = 1.0635315835416809e-08     (/mV) 
  c1h = 1.0255638036875913e-09     (/mV2) 
  d1h = -4.344940777018557e-11     (/mV3) 
  b2h = 1.0624778210633314e-08     (/mV) 
  c2h = 1.0271821295948713e-09     (/mV2) 
  d2h = -4.346406835109009e-11     (/mV3) 

  an = 0.09999128853871717     (/mV) 
  bn = -4.242674611780865     (1) 
  vhn = 49.101196573147654     (mV) 
  An = 2.6597789409230477     (/ms) 
  b1n = -0.030255147982139626     (/mV) 
  c1n = 0.0010446550877891614     (/mV2) 
  d1n = -1.0867306368149552e-05     (/mV3) 
  b2n = -0.23013571183586506     (/mV) 
  c2n = -0.0036671986209911495     (/mV2) 
  d2n = -1.5566101261108886e-05     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hInf 
  hTau 
  nInf 
  nTau 
}

STATE
{
  h
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*n*n*n*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}