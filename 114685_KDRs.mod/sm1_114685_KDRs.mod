NEURON
{
  SUFFIX KDRs 
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

  apslowd = 0.12128503195737948     (/mV) 
  bpslowd = -2.665827242859634     (1) 
  vhpslowd = -32.328348641925786     (mV) 
  Apslowd = 40.86787011206162     (/ms) 
  b1pslowd = -0.06888569832512148     (/mV) 
  c1pslowd = 0.0007388393930320863     (/mV2) 
  d1pslowd = -2.57623089932197e-06     (/mV3) 
  b2pslowd = -0.06313366254444581     (/mV) 
  c2pslowd = -0.0007602899068018277     (/mV2) 
  d2pslowd = -3.7347827583415916e-06     (/mV3) 

  apslowi = -0.02672608637195443     (/mV) 
  bpslowi = -0.016065808412289377     (1) 
  vhpslowi = 56.42914934234174     (mV) 
  Apslowi = 4489.708570114549     (/ms) 
  b1pslowi = -0.06288229487055744     (/mV) 
  c1pslowi = 0.000958426546308616     (/mV2) 
  d1pslowi = -6.317423344204479e-06     (/mV3) 
  b2pslowi = -0.04231313639232091     (/mV) 
  c2pslowi = -0.0003487456034072356     (/mV2) 
  d2pslowi = -9.552961069978118e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  pslowdInf 
  pslowdTau 
  pslowiInf 
  pslowiTau 
}

STATE
{
  pslowd
  pslowi
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*pslowd*pslowi
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  pslowd' = (pslowdInf - pslowd) / pslowdTau 
  pslowi' = (pslowiInf - pslowi) / pslowiTau 
}

INITIAL
{
  rates(v)
  pslowd = pslowdInf 
  pslowi = pslowiInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    pslowdInf = 1/(1 + exp(-apslowd*v + bpslowd)) 
    pslowdTau = Apslowd / ( exp(-(b1pslowd*(v-vhpslowd) + c1pslowd*(v-vhpslowd)^2 + d1pslowd*(v-vhpslowd)^3)) + exp((b2pslowd*(v-vhpslowd) + c2pslowd*(v-vhpslowd)^2 + d2pslowd*(v-vhpslowd)^3)) ) 

    pslowiInf = 1/(1 + exp(-apslowi*v + bpslowi)) 
    pslowiTau = Apslowi / ( exp(-(b1pslowi*(v-vhpslowi) + c1pslowi*(v-vhpslowi)^2 + d1pslowi*(v-vhpslowi)^3)) + exp((b2pslowi*(v-vhpslowi) + c2pslowi*(v-vhpslowi)^2 + d2pslowi*(v-vhpslowi)^3)) ) 


  UNITSON
}