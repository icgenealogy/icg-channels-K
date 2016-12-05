TITLE K-fast channel from Korngreen and Sakmann 2000
: M.Migliore June 2001

NEURON {
	SUFFIX iA
	USEION k READ ek WRITE ik
        RANGE gbar,ninf,linf,taul,taun
        GLOBAL tq,qq, q10,vmin, vmax,tadj,temp
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S)  = (siemens)
	(um) = (micron)

}

PARAMETER {
	v (mV)
	celsius		(degC)
	Tscale = 10	(degC)
	gbar=1.0 (S/um2)
      vhalfn=-47   (mV)
      vhalfl=-66   (mV)
      kn=29   (mV)
      kl=-10   (mV)
	qq=5
	tq=-55
	ek     (mV)
	ekkai   = -80   (mV)
	vmin = -120	(mV)
	vmax = 100	(mV)
        temp = 21       (degC)          : original temp
        q10  = 2.3
}



STATE {
	n
        l
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        linf      
        taul  (ms)
        taun   (ms)
	tadj
}

INITIAL {
	rates(v)
	n=ninf
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gbar*n^4*l*(v-ekkai)
}


DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v)
        n' = (ninf - n)/taun
        l' =  (linf - l)/taul
}

PROCEDURE rates(v (mV)) { :callable from hoc

        tadj= q10^((celsius-temp)/Tscale)
        ninf = 1/(1 + exp(-(v-vhalfn)/kn))
        linf = 1/(1 + exp(-(v-vhalfl)/kl))
        taun = (0.34+0.92*exp(-((v+71)/59)^2))
        taul = (8+49*exp(-((v+73)/23)^2))
}














