NEURON
{
  SUFFIX KAsm 
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

  ah = -0.09607501914294132     (/mV) 
  bh = 7.5713349183709955     (1) 
  vhh = -39.98971379669819     (mV) 
  Ah = 818.7656209250378     (/ms) 
  b1h = -0.0015647631362799017     (/mV) 
  c1h = -0.00016521542291544466     (/mV2) 
  d1h = 1.1216315400767086e-06     (/mV3) 
  b2h = -0.12489352781303037     (/mV) 
  c2h = -0.002145623837402677     (/mV2) 
  d2h = 2.133240301083469e-05     (/mV3) 

  am = 0.07518787297155319     (/mV) 
  bm = -1.9248051772630759     (1) 
  vhm = -37.93113674828516     (mV) 
  Am = 33.24183346425697     (/ms) 
  b1m = 0.03732507060542565     (/mV) 
  c1m = 9.756534940445149e-06     (/mV2) 
  d1m = 4.247772777941525e-08     (/mV3) 
  b2m = 0.035899776971022666     (/mV) 
  c2m = 9.669863670012e-06     (/mV2) 
  d2m = -5.121461343182813e-08     (/mV3) 
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
  mInf 
  mTau 
}

STATE
{
  h
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*m
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}