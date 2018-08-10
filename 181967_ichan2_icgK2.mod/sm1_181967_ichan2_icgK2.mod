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

  ah = -0.14595341810158313     (/mV) 
  bh = 6.940828863810799     (1) 
  vhh = -44.17324004313369     (mV) 
  Ah = 9.316662816152304     (/ms) 
  b1h = -0.12156942666927618     (/mV) 
  c1h = 0.0013434042600609722     (/mV2) 
  d1h = -4.655911769330048e-06     (/mV3) 
  b2h = -0.03512504334670375     (/mV) 
  c2h = 0.00020030211749635032     (/mV2) 
  d2h = 3.005896262044371e-07     (/mV3) 

  anf = 0.12350342059069137     (/mV) 
  bnf = -3.2748071283477125     (1) 
  vhnf = -29.53930524169872     (mV) 
  Anf = 6.0947972475613605     (/ms) 
  b1nf = -0.0989766727186109     (/mV) 
  c1nf = 0.0010385435117861962     (/mV2) 
  d1nf = -4.001710363975052e-06     (/mV3) 
  b2nf = -0.03558086678710272     (/mV) 
  c2nf = -0.00025440472964726086     (/mV2) 
  d2nf = -1.8855387355141772e-06     (/mV3) 

  am = 0.12498507946055541     (/mV) 
  bm = -3.626153323474534     (1) 
  vhm = -34.86329658157132     (mV) 
  Am = 0.2292939553968912     (/ms) 
  b1m = -0.03773332046669612     (/mV) 
  c1m = 0.00027307118759936817     (/mV2) 
  d1m = -8.110744669882462e-07     (/mV3) 
  b2m = -0.05484467321830793     (/mV) 
  c2m = -0.0008159509976556636     (/mV2) 
  d2m = -5.772887994467185e-06     (/mV3) 

  ans = 0.1235324929649318     (/mV) 
  bns = -4.757012307777989     (1) 
  vhns = -44.733358945768714     (mV) 
  Ans = 16.599179338130533     (/ms) 
  b1ns = -0.09332783597839941     (/mV) 
  c1ns = 0.0008714973773694366     (/mV2) 
  d1ns = -2.985507991339944e-06     (/mV3) 
  b2ns = -0.05075126343889554     (/mV) 
  c2ns = -0.0007156483141593805     (/mV2) 
  d2ns = -6.376485033957365e-06     (/mV3) 
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
  g = gbar*ns*ns*ns*ns
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