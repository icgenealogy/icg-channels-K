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

  ah = -0.1458419106710784     (/mV) 
  bh = 6.935000847483515     (1) 
  vhh = -22.49703730873213     (mV) 
  Ah = 0.12705732459683797     (/ms) 
  b1h = -0.025279428032321617     (/mV) 
  c1h = -0.0003517030227079615     (/mV2) 
  d1h = 3.912077719854849e-06     (/mV3) 
  b2h = -0.16598482861443015     (/mV) 
  c2h = -0.014179660867194988     (/mV2) 
  d2h = 0.0001274257931782122     (/mV3) 

  anf = 0.12312442893902273     (/mV) 
  bnf = -3.2585017512319827     (1) 
  vhnf = -38.74187413806826     (mV) 
  Anf = 0.2549288109969444     (/ms) 
  b1nf = -0.08427033334778274     (/mV) 
  c1nf = 0.0010802792932731032     (/mV2) 
  d1nf = -4.854955637429655e-06     (/mV3) 
  b2nf = -0.07278480220074394     (/mV) 
  c2nf = -0.0013577785929711956     (/mV2) 
  d2nf = -1.1153427766665097e-05     (/mV3) 

  am = 0.12460800769483199     (/mV) 
  bm = -3.609803535109155     (1) 
  vhm = -33.4789273689445     (mV) 
  Am = 0.3104330295803249     (/ms) 
  b1m = 0.08167212399548766     (/mV) 
  c1m = 0.0014470290062768404     (/mV2) 
  d1m = 1.1571200647983319e-05     (/mV3) 
  b2m = 0.10734460899468398     (/mV) 
  c2m = -0.0018037448742064968     (/mV2) 
  d2m = 1.1156276800255896e-05     (/mV3) 

  ans = 0.12327809740934449     (/mV) 
  bns = -4.74506525003381     (1) 
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