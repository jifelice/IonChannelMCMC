# Archivo con los parámetros iniciales que se van a usar para el MCMC de Tratamiento Térmico (TT)

##---- Valores de sigma, N y número de pasos para cada ciclo

#sigmasTT <- c(0.10, 0.075, 0.06, 0.05, 0.04, 0.03, 0.025, 0.02    )		# Debe tener tantos elementos como los otros vectores
#NsTT     <- c(  40,    60,   80,  100   120,  140,  140     )					# Debe tener tantos elementos como los otros vectores
#pasosTT  <- c(5000, 10000,10000, 10000, 10000, 10000, 10000)			# Debe tener tantos elementos como los otros vectores
sigmasTT=0.01
xxGxx=   0.01
for(iiGii in 2:25){
  xxGxx=xxGxx*0.80
  sigmasTT=c(sigmasTT,xxGxx)}
NsTT=60+ (0:24)*10
pasosTT= 1:25*0 + 10000
