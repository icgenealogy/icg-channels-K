NEURON
{
  SUFFIX kad 
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

  an = 0.04761897351413687     (/mV) 
  bn = -1.6380898779079844     (1) 
  vhn = -89.81459630349289     (mV) 
  An = 0.41489588593985355     (/ms) 
  b1n = -0.0022250753656668566     (/mV) 
  c1n = -4.5987921390699366e-05     (/mV2) 
  d1n = 2.636058596575081e-07     (/mV3) 
  b2n = 0.0022273509704579574     (/mV) 
  c2n = -0.00015949796987523577     (/mV2) 
  d2n = 7.407043479797671e-07     (/mV3) 

  al = -0.12195087037234278     (/mV) 
  bl = 7.073152224682837     (1) 
  vhl = -34.92831255735326     (mV) 
  Al = 9.683803698148349     (/ms) 
  b1l = 0.006274528148983956     (/mV) 
  c1l = 8.071998231259114e-05     (/mV2) 
  d1l = -4.0814460753188786e-07     (/mV3) 
  b2l = -0.006276599771485723     (/mV) 
  c2l = -0.0006639315482368562     (/mV2) 
  d2l = -8.964660267781294e-06     (/mV3) 
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