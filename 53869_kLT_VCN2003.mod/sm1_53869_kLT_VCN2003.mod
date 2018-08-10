NEURON
{
  SUFFIX kLT_VCN2003 
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

  aw = 0.08046177255517653     (/mV) 
  bw = -5.356620829714194     (1) 
  vhw = -70.10209966783464     (mV) 
  Aw = 2.4680627930291132     (/ms) 
  b1w = -0.07485762660815536     (/mV) 
  c1w = 0.0007668683133353895     (/mV2) 
  d1w = -2.4000171198097527e-06     (/mV3) 
  b2w = -0.10037814143155624     (/mV) 
  c2w = -0.004009023492916943     (/mV2) 
  d2w = -6.455584995212921e-05     (/mV3) 

  az = -0.0071927407129752045     (/mV) 
  bz = -0.31284261394443763     (1) 
  vhz = -59.30523125481516     (mV) 
  Az = 216.1505015240804     (/ms) 
  b1z = 0.12099964161099953     (/mV) 
  c1z = 0.000749946532267953     (/mV2) 
  d1z = -9.748844247904737e-06     (/mV3) 
  b2z = 0.057636661823916566     (/mV) 
  c2z = -0.0003230943942532285     (/mV2) 
  d2z = 4.939859431055287e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  wInf 
  wTau 
  zInf 
  zTau 
}

STATE
{
  w
  z
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*w*z
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  w' = (wInf - w) / wTau 
  z' = (zInf - z) / zTau 
}

INITIAL
{
  rates(v)
  w = wInf 
  z = zInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    wInf = 1/(1 + exp(-aw*v + bw)) 
    wTau = Aw / ( exp(-(b1w*(v-vhw) + c1w*(v-vhw)^2 + d1w*(v-vhw)^3)) + exp((b2w*(v-vhw) + c2w*(v-vhw)^2 + d2w*(v-vhw)^3)) ) 

    zInf = 1/(1 + exp(-az*v + bz)) 
    zTau = Az / ( exp(-(b1z*(v-vhz) + c1z*(v-vhz)^2 + d1z*(v-vhz)^3)) + exp((b2z*(v-vhz) + c2z*(v-vhz)^2 + d2z*(v-vhz)^3)) ) 


  UNITSON
}