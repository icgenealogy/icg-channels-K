NEURON
{
  SUFFIX hh3 
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

  as = -0.3332458842669122     (/mV) 
  bs = 14.662808808133764     (1) 
  vhs = -29.969152217553166     (mV) 
  As = 50.5384369088862     (/ms) 
  b1s = -0.8618792706787135     (/mV) 
  c1s = 0.0143346376499802     (/mV2) 
  d1s = -6.229960808373416e-05     (/mV3) 
  b2s = -0.00012629111057988046     (/mV) 
  c2s = -4.009789424919966e-06     (/mV2) 
  d2s = -3.367617222301758e-08     (/mV3) 

  an = 0.3331819798135625     (/mV) 
  bn = -13.326774799834855     (1) 
  vhn = -133.7065144412249     (mV) 
  An = 1.304199906704581     (/ms) 
  b1n = -0.004436762724831339     (/mV) 
  c1n = 2.4679018154833104e-05     (/mV2) 
  d1n = -4.5346292980734016e-08     (/mV3) 
  b2n = -0.08824643085179805     (/mV) 
  c2n = 0.001270321237668506     (/mV2) 
  d2n = -7.578269775964131e-06     (/mV3) 

  ah = -0.33329280692342955     (/mV) 
  bh = 14.998089255823972     (1) 
  vhh = -86.57443933502555     (mV) 
  Ah = 1.012965698248337     (/ms) 
  b1h = -0.0014896131085512807     (/mV) 
  c1h = 7.572166192034924e-05     (/mV2) 
  d1h = -3.595018575607755e-07     (/mV3) 
  b2h = 1.0035143192149093e-05     (/mV) 
  c2h = 4.344001879689698e-05     (/mV2) 
  d2h = -2.3030874432746286e-07     (/mV3) 

  am = 0.33318197981395725     (/mV) 
  bm = -13.326774799850368     (1) 
  vhm = -99.99960163103539     (mV) 
  Am = 0.12034751993818232     (/ms) 
  b1m = -0.005166398618783373     (/mV) 
  c1m = 6.020561758015037e-05     (/mV2) 
  d1m = -2.0316532740937709e-07     (/mV3) 
  b2m = 0.005164136545652506     (/mV) 
  c2m = -6.022722872078518e-05     (/mV2) 
  d2m = 2.032205814255154e-07     (/mV3) 

  an2 = 0.33318150230743404     (/mV) 
  bn2 = -13.326755365572687     (1) 
  vhn2 = -101.10499460605429     (mV) 
  An2 = 20.022663560970678     (/ms) 
  b1n2 = -3.011934389747751e-05     (/mV) 
  c1n2 = -1.6921307077775418e-05     (/mV2) 
  d1n2 = 3.723564157825848e-08     (/mV3) 
  b2n2 = -6.953014009802839e-06     (/mV) 
  c2n2 = -1.484483591071475e-05     (/mV2) 
  d2n2 = -4.824782987902349e-09     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  sInf 
  sTau 
  nInf 
  nTau 
  hInf 
  hTau 
  mInf 
  mTau 
  n2Inf 
  n2Tau 
}

STATE
{
  s
  n
  h
  m
  n2
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*s*n*n*h*m*n2
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  s' = (sInf - s) / sTau 
  n' = (nInf - n) / nTau 
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
  n2' = (n2Inf - n2) / n2Tau 
}

INITIAL
{
  rates(v)
  s = sInf 
  n = nInf 
  h = hInf 
  m = mInf 
  n2 = n2Inf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    sInf = 1/(1 + exp(-as*v + bs)) 
    sTau = As / ( exp(-(b1s*(v-vhs) + c1s*(v-vhs)^2 + d1s*(v-vhs)^3)) + exp((b2s*(v-vhs) + c2s*(v-vhs)^2 + d2s*(v-vhs)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    n2Inf = 1/(1 + exp(-an2*v + bn2)) 
    n2Tau = An2 / ( exp(-(b1n2*(v-vhn2) + c1n2*(v-vhn2)^2 + d1n2*(v-vhn2)^3)) + exp((b2n2*(v-vhn2) + c2n2*(v-vhn2)^2 + d2n2*(v-vhn2)^3)) ) 


  UNITSON
}