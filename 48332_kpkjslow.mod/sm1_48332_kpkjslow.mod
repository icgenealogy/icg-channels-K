NEURON
{
  SUFFIX kpkjslow 
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

  an = 0.05434747699136884     (/mV) 
  bn = -1.4945502125406331     (1) 
  vhn = -25.776855155217056     (mV) 
  An = 10.439946050934449     (/ms) 
  b1n = 0.0029304459973415685     (/mV) 
  c1n = -8.503422482761329e-05     (/mV2) 
  d1n = -1.4837787890849486e-07     (/mV3) 
  b2n = 0.09259423760022166     (/mV) 
  c2n = -0.0010970582110580848     (/mV2) 
  d2n = 4.053870013187504e-06     (/mV3) 
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
  g = gbar*n
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