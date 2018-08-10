NEURON
{
  SUFFIX kht 
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

  ap = 0.16665845048476746     (/mV) 
  bp = -3.8330685392689614     (1) 
  vhp = -57.36786634171382     (mV) 
  Ap = 6.318207002331458     (/ms) 
  b1p = 0.04356442617674611     (/mV) 
  c1p = 0.0002574150807540894     (/mV2) 
  d1p = -4.932783148967089e-07     (/mV3) 
  b2p = 0.03311497279892978     (/mV) 
  c2p = -0.00020652892006093443     (/mV2) 
  d2p = 4.4839641863860016e-07     (/mV3) 

  an = 0.14733584539559136     (/mV) 
  bn = -3.113045942476906     (1) 
  vhn = -45.28183732701339     (mV) 
  An = 1.5074157689495844     (/ms) 
  b1n = 0.032048349383709325     (/mV) 
  c1n = -1.959025456127755e-05     (/mV2) 
  d1n = -1.1554784735414761e-06     (/mV3) 
  b2n = 0.05105502469104417     (/mV) 
  c2n = -0.0003790828722664421     (/mV2) 
  d2n = 9.547944601033657e-07     (/mV3) 
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
  nInf 
  nTau 
}

STATE
{
  p
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*p*n*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  p' = (pInf - p) / pTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  p = pInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    pInf = 1/(1 + exp(-ap*v + bp)) 
    pTau = Ap / ( exp(-(b1p*(v-vhp) + c1p*(v-vhp)^2 + d1p*(v-vhp)^3)) + exp((b2p*(v-vhp) + c2p*(v-vhp)^2 + d2p*(v-vhp)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}