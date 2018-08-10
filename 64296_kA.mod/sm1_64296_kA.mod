NEURON
{
  SUFFIX kA 
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

  ap = 0.07692291398869601     (/mV) 
  bp = -3.2307580403552385     (1) 
  vhp = -115.83040089796961     (mV) 
  Ap = 1.4059149976938539     (/ms) 
  b1p = -0.0005622747973105617     (/mV) 
  c1p = 4.97881915540354e-06     (/mV2) 
  d1p = -1.219449199537313e-08     (/mV3) 
  b2p = -0.42593858903916065     (/mV) 
  c2p = 0.004934769608046369     (/mV2) 
  d2p = -1.5485100602436156e-05     (/mV3) 

  aq = -0.05555358675996277     (/mV) 
  bq = 6.111127376112253     (1) 
  vhq = -38.86615384173577     (mV) 
  Aq = 300.0216058319483     (/ms) 
  b1q = -0.008727207582764355     (/mV) 
  c1q = 3.8744724195383105e-05     (/mV2) 
  d1q = -6.088551211004413e-08     (/mV3) 
  b2q = -0.008726581290561655     (/mV) 
  c2q = -3.738384932530179e-05     (/mV2) 
  d2q = -5.0293708459603075e-08     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  pInf 
  pTau 
  qInf 
  qTau 
}

STATE
{
  p
  q
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*p*q
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  p' = (pInf - p) / pTau 
  q' = (qInf - q) / qTau 
}

INITIAL
{
  rates(v)
  p = pInf 
  q = qInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    pInf = 1/(1 + exp(-ap*v + bp)) 
    pTau = Ap / ( exp(-(b1p*(v-vhp) + c1p*(v-vhp)^2 + d1p*(v-vhp)^3)) + exp((b2p*(v-vhp) + c2p*(v-vhp)^2 + d2p*(v-vhp)^3)) ) 

    qInf = 1/(1 + exp(-aq*v + bq)) 
    qTau = Aq / ( exp(-(b1q*(v-vhq) + c1q*(v-vhq)^2 + d1q*(v-vhq)^3)) + exp((b2q*(v-vhq) + c2q*(v-vhq)^2 + d2q*(v-vhq)^3)) ) 


  UNITSON
}