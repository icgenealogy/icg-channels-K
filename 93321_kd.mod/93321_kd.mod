COMMENT
This file, kd.mod, implements the IKd potassium current from 
Liu et al. 1998 (Activity dependent conductances) table p.2319
Tom M Morse 20070803
ENDCOMMENT

NEURON {
	SUFFIX kd
	:NONSPECIFIC_CURRENT i
	USEION k READ ek WRITE ik
        :POINTER gbar
	RANGE gbar
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar =1.0 (S/cm2) : = 2e-6	(S/cm2) < 0, 1e9 > : this value gets overwritten by activity dependent regulation
	:Erev = -80 (mV)
}

ASSIGNED {
        ek (mV)
	ik (mA/cm2)
	v (mV)
	g (S/cm2)
	minf
	tau_m (ms)
}

STATE {	m }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m^4
	ik = g * (v - ek)
}

INITIAL {
	: assume that v has been constant for a long time
	rates(v)
	m = minf
}
DERIVATIVE states {
	rates(v)
	m' = (minf - m)/tau_m
}

FUNCTION taum(Vm (mV)) (ms) {
	UNITSOFF
	taum = 7.2-6.4/(1+exp(-(Vm+28.3)/19.2))
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	tau_m = taum(Vm)
	UNITSOFF
	minf = 1/(1+exp(-(Vm+12.3)/11.8))
	UNITSON
}
