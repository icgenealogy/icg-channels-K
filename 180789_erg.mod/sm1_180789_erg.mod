NEURON
{
  SUFFIX erg 
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

  ah = -0.049999879356045926     (/mV) 
  bh = 3.499992182187303     (1) 
  vhh = -69.68712396686963     (mV) 
  Ah = 40.75487608459169     (/ms) 
  b1h = -0.02012446174627278     (/mV) 
  c1h = 3.333628616962877e-07     (/mV2) 
  d1h = 8.349765196859151e-10     (/mV3) 
  b2h = -0.029575744678263613     (/mV) 
  c2h = 8.254016694345055e-06     (/mV2) 
  d2h = 1.0521016146118181e-07     (/mV3) 

  an = 0.20158453097178033     (/mV) 
  bn = -7.066447246486858     (1) 
  vhn = -23.705412338979723     (mV) 
  An = 7642.079481838766     (/ms) 
  b1n = -0.11999230177632876     (/mV) 
  c1n = -3.617553578193422e-07     (/mV2) 
  d1n = 7.1492330937697855e-09     (/mV3) 
  b2n = -0.05000390979289115     (/mV) 
  c2n = -7.458758768574795e-08     (/mV2) 
  d2n = -3.9117394181425765e-10     (/mV3) 
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