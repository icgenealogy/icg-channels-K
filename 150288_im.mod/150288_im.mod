: voltage-gated persistent muscarinic channel

NEURON {
	SUFFIX im
	USEION k READ ek WRITE ik
	RANGE gbar, gm
	RANGE inf, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gbar = 0.0003 (siemens/cm2) <0,1e9>
}

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2)
	inf
	tau (ms)
	gm (siemens/cm2)
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gm = gbar*n*n
	ik = gm*(v-ek)
}

INITIAL {
	rate(v)
	n = inf
}

DERIVATIVE states {
	rate(v)
	n' = (inf-n)/tau
}

FUNCTION alf(v (mV)) (/ms) {
	UNITSOFF
	alf = 0.016/exp(-(v+52.7)/23)
	UNITSON
}

FUNCTION bet(v (mV)) (/ms) {
	UNITSOFF
	bet = 0.016/exp((v+52.7)/18.8)
	UNITSON
}

PROCEDURE rate(v (mV)) {
	LOCAL sum, aa, ab
	UNITSOFF
	aa=alf(v) ab=bet(v) 
	
	sum = aa+ab
	inf = aa/sum
	tau = 1/sum
	UNITSON
}
