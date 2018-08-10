NEURON
{
  SUFFIX kv3_gp 
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

  ac = -0.1906290484793658     (/mV) 
  bc = -0.54592766682067     (1) 
  vhc = 3.0701130737121534     (mV) 
  Ac = 7.300382278362758     (/ms) 
  b1c = -0.125550891057364     (/mV) 
  c1c = -6.002292295571658e-05     (/mV2) 
  d1c = 4.954319364894917e-06     (/mV3) 
  b2c = -0.06349617221760802     (/mV) 
  c2c = 0.0001061199147609889     (/mV2) 
  d2c = 1.3604999106228759e-06     (/mV3) 

  ao = 0.19062904846976397     (/mV) 
  bo = 0.5459276667444546     (1) 
  vho = 3.232038580753739     (mV) 
  Ao = 7.26336025437143     (/ms) 
  b1o = -0.12752032135243252     (/mV) 
  c1o = 4.440077789215226e-05     (/mV2) 
  d1o = 3.094709374165125e-06     (/mV3) 
  b2o = -0.06243967216760325     (/mV) 
  c2o = 0.00012872982861899656     (/mV2) 
  d2o = 1.5072195964660658e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ek	(mV)
  ik	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  cInf 
  cTau 
  oInf 
  oTau 
}

STATE
{
  c
  o
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*o
  ik = g*(v-ek)
}

DERIVATIVE states
{
  rates(v)
  c' = (cInf - c) / cTau 
  o' = (oInf - o) / oTau 
}

INITIAL
{
  rates(v)
  c = cInf 
  o = oInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    cInf = 1/(1 + exp(-ac*v + bc)) 
    cTau = Ac / ( exp(-(b1c*(v-vhc) + c1c*(v-vhc)^2 + d1c*(v-vhc)^3)) + exp((b2c*(v-vhc) + c2c*(v-vhc)^2 + d2c*(v-vhc)^3)) ) 

    oInf = 1/(1 + exp(-ao*v + bo)) 
    oTau = Ao / ( exp(-(b1o*(v-vho) + c1o*(v-vho)^2 + d1o*(v-vho)^3)) + exp((b2o*(v-vho) + c2o*(v-vho)^2 + d2o*(v-vho)^3)) ) 


  UNITSON
}