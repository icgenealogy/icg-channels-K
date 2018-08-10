NEURON
{
  SUFFIX Ksoma 
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

  an = 0.045145471077057894     (/mV) 
  bn = -0.7291821102756127     (1) 
  vhn = 31.11613244417689     (mV) 
  An = 3.594293793522208     (/ms) 
  b1n = -0.025794123781870244     (/mV) 
  c1n = 0.00011569216204112057     (/mV2) 
  d1n = 1.6042614049822786e-07     (/mV3) 
  b2n = 0.0018375163935172682     (/mV) 
  c2n = 8.28280681415256e-05     (/mV2) 
  d2n = 2.664790596160741e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  nInf 
  nTau 
}

STATE
{
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n*n*n*n
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}