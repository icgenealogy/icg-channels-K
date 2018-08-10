NEURON
{
  SUFFIX kas 
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

  ah = -0.046443049304931185     (/mV) 
  bh = 1.557392093482853     (1) 
  vhh = 105.19924381355239     (mV) 
  Ah = 307.836299783021     (/ms) 
  b1h = -0.003682854243684531     (/mV) 
  c1h = 9.496943509175018e-05     (/mV2) 
  d1h = 4.930533113554019e-07     (/mV3) 
  b2h = 0.01730045873902411     (/mV) 
  c2h = 0.00010542476809345525     (/mV2) 
  d2h = 5.093901140425704e-07     (/mV3) 

  am = 0.0624705420175834     (/mV) 
  bm = -1.6867429573423824     (1) 
  vhm = -17.765028040302163     (mV) 
  Am = 15.41024992320288     (/ms) 
  b1m = -0.07318921478369787     (/mV) 
  c1m = -0.00031487900247672705     (/mV2) 
  d1m = 8.248483356602032e-06     (/mV3) 
  b2m = 1.6463932468908045e-05     (/mV) 
  c2m = 0.000502148230049625     (/mV2) 
  d2m = -1.8137548593387144e-06     (/mV3) 
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
  g = gbar*h*m*m
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