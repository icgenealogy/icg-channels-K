NEURON
{
  SUFFIX borgkdr 
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

  an = 0.18705029129052508     (/mV) 
  bn = -5.985599320715162     (1) 
  vhn = -31.978247823599673     (mV) 
  An = 15.466527658027973     (/ms) 
  b1n = -0.07437188316538353     (/mV) 
  c1n = -3.153651865358343e-05     (/mV2) 
  d1n = 5.233675856486221e-07     (/mV3) 
  b2n = -0.11124451925172214     (/mV) 
  c2n = 7.218107277046221e-05     (/mV2) 
  d2n = 1.821119719013095e-06     (/mV3) 

  al = -0.07482102959401694     (/mV) 
  bl = 4.564518384282611     (1) 
  vhl = 112.0547524196466     (mV) 
  Al = 945.7018645612388     (/ms) 
  b1l = 0.0020589883184219423     (/mV) 
  c1l = -5.9312218554164875e-05     (/mV2) 
  d1l = -4.905548501229527e-07     (/mV3) 
  b2l = -0.002058701088555703     (/mV) 
  c2l = -0.00017556305080957093     (/mV2) 
  d2l = -1.1826526747790789e-06     (/mV3) 
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
  lInf 
  lTau 
}

STATE
{
  n
  l
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n*l
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
  l' = (lInf - l) / lTau 
}

INITIAL
{
  rates(v)
  n = nInf 
  l = lInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    lInf = 1/(1 + exp(-al*v + bl)) 
    lTau = Al / ( exp(-(b1l*(v-vhl) + c1l*(v-vhl)^2 + d1l*(v-vhl)^3)) + exp((b2l*(v-vhl) + c2l*(v-vhl)^2 + d2l*(v-vhl)^3)) ) 


  UNITSON
}