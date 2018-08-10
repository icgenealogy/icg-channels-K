NEURON
{
  SUFFIX IKs 
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

  ah = -0.15144370005917662     (/mV) 
  bh = 10.299155034562872     (1) 
  vhh = 16.68152282970325     (mV) 
  Ah = 1048.1730359382634     (/ms) 
  b1h = 7.550487813143508e-05     (/mV) 
  c1h = 0.00011700384250083933     (/mV2) 
  d1h = -1.0323266807274552e-06     (/mV3) 
  b2h = -7.3943687187888e-05     (/mV) 
  c2h = 7.278828569041654e-05     (/mV2) 
  d2h = -5.501114922226398e-07     (/mV3) 

  am = 0.1538469635459785     (/mV) 
  bm = -5.23076238432841     (1) 
  vhm = -76.06035615261185     (mV) 
  Am = 20.012126612502552     (/ms) 
  b1m = 0.0012536251362160332     (/mV) 
  c1m = 2.978224254043667e-06     (/mV2) 
  d1m = 8.811834598295047e-09     (/mV3) 
  b2m = 0.0012057054510999164     (/mV) 
  c2m = 2.337363314208179e-06     (/mV2) 
  d2m = -3.763421234501815e-09     (/mV3) 
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