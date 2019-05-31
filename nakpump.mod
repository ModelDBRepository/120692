: Chapman JB, Johnson EA, Kootsey JM. (1983)
: Electrical and Biochemical Properties of an Enzyme Model of the Sodium Pump
: J. Membrane Biol. 74, 139-153

: note default step 5 voltage dependence

NEURON {
	SUFFIX nakpump
	RANGE atpact
	RANGE nain, naout, kin, kout, atp, p, adp
	: following for fig 9 and 12
	RANGE rf, rb
}

UNITS {
	(l) = (liter)
	(mol) = (1)
	(mmol) = (millimol)
	(uA) = (microamp)
	(mV) = (millivolt)
	F = (faraday)  (kilocoulombs)
	R = (k-mole) (joule/degC)
}

PARAMETER {
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
}

ASSIGNED {
	v (mV)
	atpact (uA/cm2)

	volin (l)
	volout (l)
	surf (cm2)
	rf[7] (uA/cm2) rb[7] (uA/cm2)
}

STATE {
	eatp (mol/cm2)
	na3eatp (mol/cm2)
	na3ep (mol/cm2)
	ep (mol/cm2)
	k2e (mol/cm2)
	k2eatp (mol/cm2)
}

INITIAL {
	volin = 1
	volout = 1
	surf = 1
	SOLVE scheme STEADYSTATE sparse
}

BREAKPOINT {
	SOLVE scheme METHOD sparse
}

KINETIC scheme {
	LOCAL vdi, vdo, a3i, a3o, a5i, a5o, x, i
	x = F/surf*(1e9)  i = 1
	a3i = exp(a3*(1 - beta)*F*v/R/T)
	a3o = exp(-a3*beta*F*v/R/T)
	a5i = exp(a5*(1 - beta)*F*v/R/T)
	a5o = exp(-a5*beta*F*v/R/T)

	COMPARTMENT volin { nain kin atp adp p }
	COMPARTMENT volout { naout kout }
	COMPARTMENT surf*(1e3) { eatp na3eatp na3ep ep k2e k2eatp }
	~ eatp + 3 nain <-> na3eatp	(f1*surf*(1e-9), b1*surf*1e0)
		rf[i] = f_flux*x  rb[i] = b_flux*x  i = i+1
	~ na3eatp <-> na3ep + adp	(f2*surf*(1e0), b2*surf*(1e-3))
		rf[i] = f_flux*x  rb[i] = b_flux*x  i = i+1
	~ na3ep <-> ep + 3 naout	(a3i*f3*surf*(1e0), a3o*b3*surf*(1e-9))
		rf[i] = f_flux*x  rb[i] = b_flux*x  i = i+1
	~ ep + 2 kout <-> k2e + p	(f4*surf*(1e-6), b4*surf*(1e-3))
		rf[i] = f_flux*x  rb[i] = b_flux*x  i = i+1
	~ k2e + atp <-> k2eatp		(a5i*f5*surf*(1e-3), a5o*b5*surf*(1e0))
		rf[i] = f_flux*x  rb[i] = b_flux*x  i = i+1
	~ k2eatp <-> eatp + 2 kin	(f6*surf*(1e0), b6*surf*(1e-6))
		rf[i] = f_flux*x  rb[i] = b_flux*x  i = i+1
	atpact = (f_flux - b_flux)*x
	CONSERVE eatp+na3eatp+na3ep+ep+k2e+k2eatp = totalpump*surf*(1e3)
}
