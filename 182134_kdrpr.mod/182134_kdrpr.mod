COMMENT

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
//
// Copyright 2007, The University Of Pennsylvania
// 	School of Engineering & Applied Science.
//   All rights reserved.
//   For research use only; commercial use prohibited.
//   Distribution without permission of Maciej T. Lazarewicz not permitted.
//   mlazarew@seas.upenn.edu
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ENDCOMMENT


NEURON {

	SUFFIX kdrpr
	USEION k READ ek WRITE ik
	RANGE gkdr, ik
	GLOBAL ek
}
	
UNITS {

	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {

    gkdr =  15 (mS/cm2)
    ek   = -75 (mV)
}
    
ASSIGNED {

    v       (mV)
    ik      (mA/cm2)
    ninf    (1)
    taun    (ms)
}

STATE { n }

INITIAL { 

    rates(v)
    n  = ninf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	
	ik = 1 * gkdr * n * (v-ek)
}


DERIVATIVE states { 

    rates(v)
    n' = (ninf-n)/taun
}

PROCEDURE rates(v(mV)) { LOCAL a, b

    a = fun3(v,  -24.9, -0.016,   -5)
    b = fun1(v,  -40,    0.25,   -40)
    
    ninf = a/(a+b)
    taun = 1.0/(a+b)
}

INCLUDE "custom_code/inc_files/182134_aux_fun.inc"
