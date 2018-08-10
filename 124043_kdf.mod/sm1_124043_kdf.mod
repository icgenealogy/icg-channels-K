NEURON
{
  SUFFIX kdf 
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

  ah = -0.09082417364326006     (/mV) 
  bh = 4.905235145538766     (1) 
  vhh = -23.4705244709602     (mV) 
  Ah = 518.0687497665291     (/ms) 
  b1h = -0.004361520350569842     (/mV) 
  c1h = -0.0002665233581788834     (/mV2) 
  d1h = 1.6447319204587027e-06     (/mV3) 
  b2h = 0.029434649302159478     (/mV) 
  c2h = -9.53659491000034e-05     (/mV2) 
  d2h = -3.074285271041627e-06     (/mV3) 

  an = 0.07140199321292454     (/mV) 
  bn = -1.7136433096959653     (1) 
  vhn = -60.596054418989674     (mV) 
  An = 23.492551228131568     (/ms) 
  b1n = 0.09346119397684943     (/mV) 
  c1n = 0.0027564839171022857     (/mV2) 
  d1n = 3.751506345910895e-05     (/mV3) 
  b2n = 0.033223271259755785     (/mV) 
  c2n = -0.00010578115416437985     (/mV2) 
  d2n = 2.4665619005153016e-07     (/mV3) 
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