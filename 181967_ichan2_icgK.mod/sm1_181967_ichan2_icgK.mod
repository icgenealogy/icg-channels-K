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

  ah = -0.14591438859926284     (/mV) 
  bh = 6.938760673501898     (1) 
  vhh = -44.14442847623243     (mV) 
  Ah = 9.293363568579673     (/ms) 
  b1h = -0.1214919636515275     (/mV) 
  c1h = 0.0013425957236243165     (/mV2) 
  d1h = -4.652893873073303e-06     (/mV3) 
  b2h = -0.035009961139243345     (/mV) 
  c2h = 0.00020202859252930203     (/mV2) 
  d2h = 2.948862958562447e-07     (/mV3) 

  anf = 0.12349521067610675     (/mV) 
  bnf = -3.2733077913276163     (1) 
  vhnf = -29.437882655671864     (mV) 
  Anf = 6.078285796361462     (/ms) 
  b1nf = -0.09898203015091585     (/mV) 
  c1nf = 0.0010387528727149365     (/mV2) 
  d1nf = -4.0020768079006396e-06     (/mV3) 
  b2nf = -0.035103288379368235     (/mV) 
  c2nf = -0.00023996411169107277     (/mV2) 
  d2nf = -1.756655181380855e-06     (/mV3) 

  am = 0.12495991334188407     (/mV) 
  bm = -3.624205941512621     (1) 
  vhm = -34.86329658157132     (mV) 
  Am = 0.2292939553968912     (/ms) 
  b1m = -0.03773332046669612     (/mV) 
  c1m = 0.00027307118759936817     (/mV2) 
  d1m = -8.110744669882462e-07     (/mV3) 
  b2m = -0.05484467321830793     (/mV) 
  c2m = -0.0008159509976556636     (/mV2) 
  d2m = -5.772887994467185e-06     (/mV3) 

  ans = 0.12347622970773128     (/mV) 
  bns = -4.754177317319358     (1) 
  vhns = -44.68962935512492     (mV) 
  Ans = 16.577114677562502     (/ms) 
  b1ns = -0.09324507778955009     (/mV) 
  c1ns = 0.0008697308323231864     (/mV2) 
  d1ns = -2.9762284299008882e-06     (/mV3) 
  b2ns = -0.05055405968582514     (/mV) 
  c2ns = -0.0007082322759409731     (/mV2) 
  d2ns = -6.3061347526833094e-06     (/mV3) 
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
  g = gbar*nf*nf*nf*nf
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