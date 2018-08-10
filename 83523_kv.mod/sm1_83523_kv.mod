NEURON
{
  SUFFIX kv 
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

  ah = -0.09083884456243303     (/mV) 
  bh = 1.635898296638295     (1) 
  vhh = 7.861141699076644     (mV) 
  Ah = 600.022696644455     (/ms) 
  b1h = 0.0001698493241223289     (/mV) 
  c1h = -2.5271009370698045e-07     (/mV2) 
  d1h = -2.1125528515947126e-08     (/mV3) 
  b2h = 0.00016948027169380744     (/mV) 
  c2h = -2.5199851811442997e-07     (/mV2) 
  d2h = -2.1077620199626422e-08     (/mV3) 

  an = 0.09090892126973449     (/mV) 
  bn = -1.3636272361247566     (1) 
  vhn = -25.212704378460764     (mV) 
  An = 5.883717929666873     (/ms) 
  b1n = -0.06150099695185009     (/mV) 
  c1n = 0.0005278606664745815     (/mV2) 
  d1n = -1.658953281725473e-06     (/mV3) 
  b2n = -0.07168133951515274     (/mV) 
  c2n = -0.0008113734087192367     (/mV2) 
  d2n = -3.732516955176238e-06     (/mV3) 
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
  g = gbar*h*n
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