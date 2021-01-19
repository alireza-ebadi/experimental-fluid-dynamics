### calculate the air properties for atmospheric pressure
def airProps(T):
	'''
	calculate the air properties for atmospheric pressure

	input:
    T (C)

	outputs:
    rho (kg/m3)
	mu (Pa*s = Kg/m/s)
	nu (m2/s)
	'''

	R = 286.9;          #J/kg*K
	T = T+273.15;       #K
	P = 101.3e3;         #Pa
	rho = P/(R*T);      #kg/m^3
	muo = 18.27e-6;     #Pa*s
	To = 291.15;        #K
	C = 120;
	mu = muo*((To+C)/(T+C))*((T/To)**1.5);   #Pa*s
	nu = mu/rho;        #m^2/s

	return rho, mu, nu
