# Programa con nucleo para diagonalizar la matriz en el proceso de los ajustes del modelo
# En este caso no es necesario ingresar el pHi (los valores a 0 V son justamente los que se ajustan)
#
# Inputs: voltajes, T, p0, tiempos
#
# Outputs: p_open_totales (vector con la variacion de estados a lo largo del tiempo para todos los voltajes)


calcularptot <- function(voltajes, q, T, p0, tiempos){

	nEstados <- length(p0)

	p_open_totales <- matrix(0, nrow=length(tiempos), ncol=length(voltajes)+1)	# Una fila será para los t también

	t_tot <- c((1:length(tiempos))*0)


	ivoltajes <- 1

	for (V in voltajes){
	
		W <- calcularW(q, V, T, nEstados)		# Calcula W para ese V

		auto_vec_val <- eigen(W) 

		auto_vectores <- auto_vec_val$vectors
		auto_valores <- auto_vec_val$values

		cj <-solve(auto_vectores, p0)

		p_tot <- matrix(0, nrow=length(tiempos), ncol=nEstados)

		p_tot[1,] <- p0

		t <- 0

		for (it in 1:length(tiempos)){				# GUARDA LOS 4 ESTADOS... ¿SOLO EL O?
			t <- tiempos[it]
			p <- c((1:nEstados)*0)
			for (i in 1:nEstados){
				p <- p + cj[i]*exp(auto_valores[i]*t)*auto_vectores[,i]
				}
			p_tot[it,] <- p
			}


		p_open_totales[,ivoltajes+1] <- p_tot[,4] # Para guardar los valores a todos los V

		ivoltajes <- ivoltajes + 1
		}
	
	p_open_totales[,1] <- tiempos # Para guardar los valores de tiempos

	return(p_open_totales)
	}