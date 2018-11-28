TITLE K-A channel from Klee Ficker and Heinemann
: modified by Brannon and Yiota Poirazi (poirazi@LNC.usc.edu) 
: to account for Hoffman et al 1997 distal region kinetics
: used only in locations > 100 microns from the soma


UNITS {
	  (mA) = (milliamp)
	  (mV) = (millivolt)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

PARAMETER { :parameters that can be entered when function is called in cell-setup   
	  v               (mV)
    ek = -80        (mV) :K reversal potential  (reset in cell-setup.hoc)
	  celsius = 24	  (degC)
    :	gkabar = 0.008  (mho/cm2)  :suggested conductance value
	  gkabar = 1.0      (mho/cm2)  :initialized conductance
    vhalfn = -1     (mV)       :activation half-potential
    vhalfl = -56    (mV)       :inactivation half-potential
    a0n = 0.1       (/ms)      :parameters used
    zetan = -1.8    (1)        :in calculation of
    zetal = 3       (1)        :steady state values
    gmn   = 0.39    (1)        :and time constants
    gml   = 1       (1)
	  lmin  = 2       (ms)
	  nmin  = 0.1     (ms)
	  pw    = -1      (1)
	  tq    = -40     (mV)
	  qq    = 5       (mV)
	  q10   = 5                  :temperature sensitivity
}

NEURON {
	  SUFFIX kad
	  USEION k READ ek WRITE ik
    RANGE gkabar,gka, gmax
    GLOBAL ninf,linf,taul,taun,lmin
}

STATE {       :the unknown parameters to be solved in the DEs 
	  n l
}

ASSIGNED {       :parameters needed to solve DE
	  ik   (mA/cm2)
    ninf
    linf      
    taul (ms)
    taun (ms)
    gka  (mho/cm2)
    gmax (mho/cm2)
}

INITIAL {		:initialize the following parameter using rates()
	  rates(v)
	  n = ninf
	  l = linf
	  gka = gkabar*n*l
	  ik = gka*(v-ek)	
    gmax = gka
}

BREAKPOINT {
	  SOLVE states  METHOD cnexp
	  gka = gkabar*n*l
	  ik = gka*(v-ek)
    if (gka > gmax) {
        gmax = gka
    }
}

FUNCTION alpn(v(mV)) { LOCAL zeta
    zeta = zetan+pw/(1+exp((v-tq)/qq))
    alpn = exp((1.e-3)*zeta*(v-vhalfn)*FARADAY/(R*(273.16(degC)+celsius))) 
}

FUNCTION betn(v(mV)) { LOCAL zeta
    zeta = zetan+pw/(1+exp((v-tq)/qq))
    betn = exp((1.e-3)*zeta*gmn*(v-vhalfn)*FARADAY/(R*(273.16(degC)+celsius))) 
}

FUNCTION alpl(v(mV)) {
    alpl = exp((1.e-3)*zetal*(v-vhalfl)*FARADAY/(R*(273.16(degC)+celsius))) 
}

FUNCTION betl(v(mV)) {
    betl = exp((1.e-3)*zetal*gml*(v-vhalfl)*FARADAY/(R*(273.16(degC)+celsius))) 
}

:if state_borgka is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
DERIVATIVE states {     : exact when v held constant; integrates over dt step
    rates(v)
    n' = (ninf - n)/taun
    l' = (linf - l)/taul
}

PROCEDURE rates(v (mV)) {                  :callable from hoc
    LOCAL a,qt
    TABLE ninf, taun, linf, taul  DEPEND celsius, vhalfn, vhalfl, a0n, zetan,  zetal, gmn, gml, lmin, nmin,	pw,	tq,	qq,	q10 FROM -100 TO 100 WITH 200
    qt = q10^((celsius-24(degC))/10(degC))         : temprature adjastment factor
    a = alpn(v)
    ninf = 1/(1 + a)                   : activation variable steady state value
    taun = betn(v)/(qt*a0n*(1+a))      : activation variable time constant
	  if (taun<nmin) {taun=nmin}         : time constant not allowed to be less than nmin
    a = alpl(v)
    linf = 1/(1+ a)                    : inactivation variable steady state value
	  taul = 0.26(ms/mV)*(v+50(mV))                 : inactivation variable time constant
	  if (taul<lmin) {taul=lmin}         : time constant not allowed to be less than lmin
}
