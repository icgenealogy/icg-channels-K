NEURON
{
  SUFFIX Kv4 
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

  ah = -0.10403299672353487     (/mV) 
  bh = 7.186494881244829     (1) 
  vhh = -25.581047841613195     (mV) 
  Ah = 9.314171087143317     (/ms) 
  b1h = -0.008339670043422683     (/mV) 
  c1h = 0.0004911769238645753     (/mV2) 
  d1h = -2.9365395886906493e-06     (/mV3) 
  b2h = 0.008340021978442727     (/mV) 
  c2h = 1.469768946794512e-05     (/mV2) 
  d2h = -3.623101092012973e-07     (/mV3) 

  an = 0.057713128600957474     (/mV) 
  bn = -3.289646993082159     (1) 
  vhn = -77.25622230966673     (mV) 
  An = 1.097003501203185     (/ms) 
  b1n = -0.011018311954205438     (/mV) 
  c1n = -0.00018600527799436344     (/mV2) 
  d1n = 6.313903633361424e-07     (/mV3) 
  b2n = -0.035952349240001534     (/mV) 
  c2n = -0.00014295988338517204     (/mV2) 
  d2n = -3.500880259688354e-06     (/mV3) 
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
  nInf 
  nTau 
}

STATE
{
  h
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}