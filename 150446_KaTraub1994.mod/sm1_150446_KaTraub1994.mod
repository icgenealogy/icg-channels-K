NEURON
{
  SUFFIX Ka 
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

  ah = -0.2470102984685315     (/mV) 
  bh = 18.093030592818714     (1) 
  vhh = -77.238566509004     (mV) 
  Ah = 965.2073999143157     (/ms) 
  b1h = -0.1697587168057256     (/mV) 
  c1h = 0.0018914764638698466     (/mV2) 
  d1h = -6.16018900200648e-06     (/mV3) 
  b2h = -0.13846819920901773     (/mV) 
  c2h = -0.005829378707951113     (/mV2) 
  d2h = -0.00013075298143224283     (/mV3) 

  am = 0.07998962680527627     (/mV) 
  bm = -3.205550280159502     (1) 
  vhm = -37.34736264634554     (mV) 
  Am = 2.916031020753281     (/ms) 
  b1m = 0.02474408290916359     (/mV) 
  c1m = -2.362776007463877e-05     (/mV2) 
  d1m = -1.3351554048015117e-06     (/mV3) 
  b2m = 0.0314958347751929     (/mV) 
  c2m = -0.00012241460044662392     (/mV2) 
  d2m = -2.2284234508025972e-07     (/mV3) 
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