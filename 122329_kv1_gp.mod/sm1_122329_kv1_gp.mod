NEURON
{
  SUFFIX kv1_gp 
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

  ah = -0.08050106211831172     (/mV) 
  bh = 4.013607732460649     (1) 
  vhh = -76.03067305835305     (mV) 
  Ah = 395.71217749522623     (/ms) 
  b1h = 0.0020335616762175877     (/mV) 
  c1h = 0.00014173171782312654     (/mV2) 
  d1h = -8.462472273060407e-07     (/mV3) 
  b2h = -0.002033506704988678     (/mV) 
  c2h = 0.00013020797961276292     (/mV2) 
  d2h = -5.745955501829356e-07     (/mV3) 

  am = 0.08264133134791192     (/mV) 
  bm = -1.9832306576132965     (1) 
  vhm = -17.543780439032314     (mV) 
  Am = 26.48172439104839     (/ms) 
  b1m = -0.07367709706721659     (/mV) 
  c1m = -0.0003118875679351066     (/mV2) 
  d1m = 8.33564243732411e-06     (/mV3) 
  b2m = 0.0005375460290349272     (/mV) 
  c2m = 0.0005024089775491553     (/mV2) 
  d2m = -1.8557683480268048e-06     (/mV3) 
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