: Chapman JB, Johnson EA, Kootsey JM. (1983)
: Electrical and Biochemical Properties of an Enzyme Model of the Sodium Pump
: J. Membrane Biol. 74, 139-153

: note default step 5 voltage dependence

: extended by Michael Hines as a component of larger models.
: I.e. modifies nai, ki, contributes to ina, ik, and consumes atp
: for investigation of isolated pump, allow clamping of
: nai, ki, atp (note p and adp are constant here)
: initialize to steady state pump with nai, ki, atp clamped.

NEURON {
	SUFFIX nakpump
	USEION na READ ina WRITE nai, nao, ina
	USEION k READ ik WRITE ki, ko, ik
	RANGE inapump, ikpump
	RANGE atpact
	RANGE nain, naout, kin, kout, atp, p, adp
	: following for fig 9 and 12
	RANGE rf, rb
}

UNITS {
	(l) = (liter)
	(mol) = (1)
	(mmol) = (millimol)
	(mM) = (mmol/l)
	(uA) = (microamp)
	(mA) = (milliamp)
	(mV) = (millivolt)
	(um) = (micron)
	F = (faraday)  (kilocoulombs)
	R = (k-mole) (joule/degC)
	PI = (pi) (1)
}

PARAMETER {
	nasrcrate = 0 (/ms)
	ksrcrate = 0 (/ms)
	atpsrcrate = 1e9 (/ms)

	totalpump = 1.25e-13 (mol/cm2)

	nain = 9.6 (mmol/l)
	naout = 140 (mmol/l)
	kin = 150.4 (mmol/l)
	kout = 5.4 (mmol/l)
	atp = 4.99 (mmol/l)
	p = 4.95 (mmol/l)
	adp = 0.06 (mmol/l)
	T = 310 (K)

	f1 = 2.5e11 (l3/mol3-s)
	b1 = 1e5 (/s)
	f2 = 1e4 (/s)
	b2 = 1e5 (l/mol-s)
	f3 = 172 (/s)
	b3 = 1.72e4 (l3/mol3-s)
	f4 = 1.5e7 (l2/mol2-s)
	b4 = 2e5 (l/mol-s)
	f5 = 2e6 (l/mol-s)
	b5 = 30 (/s)
	f6 = 1.15e4 (/s)
	b6 = 6e8 (l2/mol2-s)

	beta = .5
	a3 = 0 (1)
	a5 = 1 (1)

	iter=1
}

ASSIGNED {
	diam (um)
	v (mV)
	ina (mA/cm2)
	ik (mA/cm2)
	inapump (mA/cm2)
	ikpump (mA/cm2)
	inapumplast (mA/cm2)
	ikpumplast (mA/cm2)
	atpact (uA/cm2)

	rf[7] (uA/cm2) rb[7] (uA/cm2)

}

STATE {
	eatp (mol/cm2)
	na3eatp (mol/cm2)
	na3ep (mol/cm2)
	ep (mol/cm2)
	k2e (mol/cm2)
	k2eatp (mol/cm2)

	nai (mM)
	ki (mM)
	nao (mM)
	ko (mM)
	atps (mmol/l)
}

LOCAL volin, volout, surf

INITIAL {
	LOCAL a1, a2, a3, a4, a5
	volin = PI*diam*diam/4 : cross section area
	volout = 1 (um2)
	surf = PI*diam*(1e7) : circumference
	nai = nain
	ki = kin
	nao = naout
	ko = kout
	atps = atp

	: clamp to nain, kin, atp
	a1 = nasrcrate  a2 = ksrcrate  a3 = atpsrcrate
	nasrcrate=1e9 ksrcrate=1e9 atpsrcrate=1e9

	ina = 0
	ik = 0
	inapump = 0
	ikpump = 0
	inapumplast = 0
	ikpumplast = 0
	FROM i = 1 TO iter {
	SOLVE scheme STEADYSTATE sparse
	}

	: unclamp
	nasrcrate=a1 ksrcrate=a2 atpsrcrate=a3
}

BREAKPOINT {
	SOLVE scheme METHOD sparse
	inapumplast = inapump
	ikpumplast = ikpump
	inapump = 3*atpact*(1e-3)
	ikpump = -2*atpact*(1e-3)
	ina = inapump
	ik = ikpump
}

KINETIC scheme {
	LOCAL vdi, vdo, a3i, a3o, a5i, a5o, x, i
	x = F/surf*(1e9)  i = 1
	a3i = exp(a3*(1 - beta)*F*v/R/T)
	a3o = exp(-a3*beta*F*v/R/T)
	a5i = exp(a5*(1 - beta)*F*v/R/T)
	a5o = exp(-a5*beta*F*v/R/T)

	COMPARTMENT volin { nai ki atps adp p }
	COMPARTMENT volout { nao ko }
	COMPARTMENT surf*(1e3) { eatp na3eatp na3ep ep k2e k2eatp }
	~ eatp + 3 nai <-> na3eatp	(f1*surf*(1e-9), b1*surf*(1e0))
		rf[i] = f_flux*x  rb[i] = b_flux*x  i = i+1
	~ na3eatp <-> na3ep + adp	(f2*surf*(1e0), b2*surf*(1e-3))
		rf[i] = f_flux*x  rb[i] = b_flux*x  i = i+1
	~ na3ep <-> ep + 3 naout	(a3i*f3*surf*(1e0), a3o*b3*surf*(1e-9))
		rf[i] = f_flux*x  rb[i] = b_flux*x  i = i+1
	~ ep + 2 kout <-> k2e + p	(f4*surf*(1e-6), b4*surf*(1e-3))
		rf[i] = f_flux*x  rb[i] = b_flux*x  i = i+1
	~ k2e + atps <-> k2eatp		(a5i*f5*surf*(1e-3), a5o*b5*surf*(1e0))
		rf[i] = f_flux*x  rb[i] = b_flux*x  i = i+1
	~ k2eatp <-> eatp + 2 ki	(f6*surf*(1e0), b6*surf*(1e-6))
		rf[i] = f_flux*x  rb[i] = b_flux*x  i = i+1
	atpact = (f_flux - b_flux)*x
	CONSERVE eatp+na3eatp+na3ep+ep+k2e+k2eatp = totalpump*surf*(1e3)

	: sources
	COMPARTMENT volin { nain kin atp }
	~ nain <-> nai (nasrcrate*volin, nasrcrate*volin)
	~ kin <-> ki (ksrcrate*volin, ksrcrate*volin)
	~ atp <-> atps (atpsrcrate*volin, atpsrcrate*volin)
	~ nai << (-(ina - inapumplast)/x*(1e3))
	~ ki << (-(ik - ikpumplast)/x*(1e3))
	COMPARTMENT volout { naout kout }
	~ naout <-> nao (1e9(/ms)*volout, 1e9(/ms)*volout)
	~ kout <-> ko (1e9(/ms)*volout, 1e9(/ms)*volout)
}
