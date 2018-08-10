NEURON
{
  SUFFIX iapnew 
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

  am = 0.3331592784389013     (/mV) 
  bm = -13.325622223023927     (1) 
  vhm = -99.99947503065728     (mV) 
  Am = 0.12258643415036646     (/ms) 
  b1m = -0.0066322389660997925     (/mV) 
  c1m = 7.062837542469439e-05     (/mV2) 
  d1m = -2.206123842229438e-07     (/mV3) 
  b2m = 0.006631614174210335     (/mV) 
  c2m = -7.063493718775373e-05     (/mV2) 
  d2m = 2.20610417307527e-07     (/mV3) 

  an = 0.3331592784468021     (/mV) 
  bn = -13.325622223227008     (1) 
  vhn = -105.19227502935148     (mV) 
  An = 2.0168215613139564     (/ms) 
  b1n = 1.279932283473638     (/mV) 
  c1n = -0.015546595440884358     (/mV2) 
  d1n = 4.8460603333885056e-05     (/mV3) 
  b2n = 0.0002911750522675683     (/mV) 
  c2n = -2.9018209583029925e-06     (/mV2) 
  d2n = 7.661068091960044e-09     (/mV3) 

  ah = -0.333266549364841     (/mV) 
  bh = 14.996827239750107     (1) 
  vhh = -155.28256420910307     (mV) 
  Ah = 0.9556489023127754     (/ms) 
  b1h = 0.011847449005285494     (/mV) 
  c1h = -0.00010061825339903402     (/mV2) 
  d1h = 5.356997351025545e-07     (/mV3) 
  b2h = 0.0037159221198170395     (/mV) 
  c2h = 9.680398752405582e-09     (/mV2) 
  d2h = -1.951413791850215e-08     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
  nInf 
  nTau 
  hInf 
  hTau 
}

STATE
{
  m
  n
  h
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*m*n*n*h
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
  n' = (nInf - n) / nTau 
  h' = (hInf - h) / hTau 
}

INITIAL
{
  rates(v)
  m = mInf 
  n = nInf 
  h = hInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 


  UNITSON
}