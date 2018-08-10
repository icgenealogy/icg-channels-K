NEURON
{
  SUFFIX IKv 
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

  ah = -0.1513735352560179     (/mV) 
  bh = 13.079256353670067     (1) 
  vhh = -81.82598085222564     (mV) 
  Ah = 0.022030113036703915     (/ms) 
  b1h = 0.24760977455472014     (/mV) 
  c1h = -0.0029826829165159156     (/mV2) 
  d1h = 9.890428889703753e-06     (/mV3) 
  b2h = 0.0037435761766186817     (/mV) 
  c2h = -4.2521535672715313e-05     (/mV2) 
  d2h = 1.373486573131778e-07     (/mV3) 

  am = 0.10079422022729094     (/mV) 
  bm = -5.221376713895346     (1) 
  vhm = -49.43753306359081     (mV) 
  Am = 0.021199262239467272     (/ms) 
  b1m = -0.00020061755475471175     (/mV) 
  c1m = -0.0002235184170512805     (/mV2) 
  d1m = -1.684009317333767e-05     (/mV3) 
  b2m = -0.289491677373407     (/mV) 
  c2m = -0.07596012822318969     (/mV2) 
  d2m = -0.004713030362539296     (/mV3) 
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
  g = gbar*h*m*m*m
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