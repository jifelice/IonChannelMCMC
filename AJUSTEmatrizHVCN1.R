# Calcula la matriz W para el conjunto de parámetros que se ingresan en el vector q
#
# Inputs: q (valores de las constantes a 0 V), V, T, nEstados
#
# Outputs: W_HVCN1 (matriz de transiciones entre los estados)


calcularW <- function(q, V, T, nEstados){		

	k <- 1.3806488e-23			# [J/K] Constante de Boltzmann (https://physics.nist.gov/cuu/Constants/index.html)
	q_electron <- 1.602e-19			# [Coulomb] Carga del electrón
	

	# Asigna la posicion en el vector q de las constantes a 0 V
	
	alpha_0_1 <- q[1]				# [mseg-1]
	alpha_0_2 <- q[2]				# [mseg-1]
	alpha_0_3 <- q[3]				# [mseg-1]
	beta_0_1 <- q[4]				# [mseg-1]
	beta_0_2 <- q[5]				# [mseg-1]
	beta_0_3 <- q[6]				# [mseg-1]
	zeta_alpha_1 <- q[7]			# [número electrones]
	zeta_alpha_2 <- q[8]			# [número electrones]
	zeta_alpha_3 <- q[9]			# [número electrones]
	zeta_beta_1 <- q[10]			# [número electrones]
	zeta_beta_2 <- q[11]			# [número electrones]
	zeta_beta_3 <- q[12]			# [número electrones]


	# Calcula las constantes (unidades de alphas y betas, mseg-1)
	
	alpha_1 <- alpha_0_1*exp((zeta_alpha_1*V*q_electron)/(k*T*1e3))
	alpha_2 <- alpha_0_2*exp((zeta_alpha_2*V*q_electron)/(k*T*1e3))
	alpha_3 <- alpha_0_3*exp((zeta_alpha_3*V*q_electron)/(k*T*1e3))
	beta_1 <- beta_0_1*exp(-(zeta_beta_1*V*q_electron)/(k*T*1e3))
	beta_2 <- beta_0_2*exp(-(zeta_beta_2*V*q_electron)/(k*T*1e3))
	beta_3 <- beta_0_3*exp(-(zeta_beta_3*V*q_electron)/(k*T*1e3))


	# Calcula la matriz a partir de las constantes calculadas

	W_HVCN1 <- matrix(c(-alpha_1, beta_1, 0, 0,
			alpha_1, -(beta_1+alpha_2), beta_2, 0,
			0, alpha_2, -(beta_2+alpha_3), beta_3, 
			0, 0, alpha_3, -beta_3), 
			nrow=nEstados, byrow=T)

	return(W_HVCN1)
	
	}