NEURON
{
  SUFFIX ichan2 
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

  ah = -0.14590168132312634     (/mV) 
  bh = 6.938120069688761     (1) 
  vhh = -59.290257501934995     (mV) 
  Ah = 0.35429936297395276     (/ms) 
  b1h = 0.1193882270330621     (/mV) 
  c1h = 0.003100024204365057     (/mV2) 
  d1h = 4.127234694333667e-05     (/mV3) 
  b2h = 0.05788499803320327     (/mV) 
  c2h = -0.0003673832225446087     (/mV2) 
  d2h = 7.308381634901477e-07     (/mV3) 

  anf = 0.12323817670243695     (/mV) 
  bnf = -3.2643284396984704     (1) 
  vhnf = -25.43918915703582     (mV) 
  Anf = 0.18398522406760032     (/ms) 
  b1nf = -0.10653614066809221     (/mV) 
  c1nf = 0.001976160275610243     (/mV2) 
  d1nf = -1.5247984159515117e-05     (/mV3) 
  b2nf = -0.017125092194956566     (/mV) 
  c2nf = 0.0001663843245489111     (/mV2) 
  d2nf = 1.5887675537941744e-06     (/mV3) 

  am = 0.12474042249422158     (/mV) 
  bm = -3.616151102378712     (1) 
  vhm = -27.860546661764225     (mV) 
  Am = 0.18591134677875681     (/ms) 
  b1m = 0.11271258954122713     (/mV) 
  c1m = 0.002265020928389483     (/mV2) 
  d1m = 1.8211282345657036e-05     (/mV3) 
  b2m = 0.12306060234245304     (/mV) 
  c2m = -0.0029523997026349967     (/mV2) 
  d2m = 3.099905918541683e-05     (/mV3) 

  ans = 0.12340020722916381     (/mV) 
  bns = -4.750866236548475     (1) 
  vhns = -45.325474447022984     (mV) 
  Ans = 0.5839592696636453     (/ms) 
  b1ns = 0.053320549702270195     (/mV) 
  c1ns = 0.0008200031204979145     (/mV2) 
  d1ns = 7.555821782871203e-06     (/mV3) 
  b2ns = 0.09288995599331366     (/mV) 
  c2ns = -0.0009188444084740952     (/mV2) 
  d2ns = 3.183772904054908e-06     (/mV3) 
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
  nfInf 
  nfTau 
  mInf 
  mTau 
  nsInf 
  nsTau 
}

STATE
{
  h
  nf
  m
  ns
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*nf*nf*nf*nf*m*m*m*ns*ns*ns*ns
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  nf' = (nfInf - nf) / nfTau 
  m' = (mInf - m) / mTau 
  ns' = (nsInf - ns) / nsTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  nf = nfInf 
  m = mInf 
  ns = nsInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nfInf = 1/(1 + exp(-anf*v + bnf)) 
    nfTau = Anf / ( exp(-(b1nf*(v-vhnf) + c1nf*(v-vhnf)^2 + d1nf*(v-vhnf)^3)) + exp((b2nf*(v-vhnf) + c2nf*(v-vhnf)^2 + d2nf*(v-vhnf)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    nsInf = 1/(1 + exp(-ans*v + bns)) 
    nsTau = Ans / ( exp(-(b1ns*(v-vhns) + c1ns*(v-vhns)^2 + d1ns*(v-vhns)^3)) + exp((b2ns*(v-vhns) + c2ns*(v-vhns)^2 + d2ns*(v-vhns)^3)) ) 


  UNITSON
}