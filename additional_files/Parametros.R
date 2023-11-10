# Archivo con los parámetros iniciales que se van a usar


##--- Parámetros generales

nPasosBase <-10000			 	# Número mínimo de pasos que se ejecuta para un sigma y un N (40k fijo, 25k variable)
fPasos <- 0					# Pasos ctes (=0) o variables (=1).
nMaxCiclos <- 15			# Número máximo de ciclos de nPasos que ejecuta
cargaMin <- 5				# Mínimo valor que puede tomar la suma de las cargas (los 6 Zs)
cargaMax <- 6				# Máximo valor que puede tomar la suma de las cargas (los 6 Zs)
told2 <- 0.2				# Tolerancia de distancia disminuida
d2Fin <- 1.e-7     			# Diferencia a la que termina el programa (break)
nPromd2 <- 1500				# numero de pasos en que promediara el valor de d2
Nini <- 20					# Valor inicial que toma el N (el mismo para todos los parámetros)
Nmax <- 40  					# Valor máximo al que se deja llegar a N
maxCiclosNSigma <- 5			# Número máximo de ciclos en que mantiene sigma y N sin variarlos
factorSigma <- 1.5			# Factor por el que se divide a sigma cada vez que cambia
factorN <- 2				# Factor por el que se multiplica a N cada vez que cambia
semilla <- XXXsem			# Establece el valor de la semilla
set.seed(semilla)				# Semilla para que los números aleatorios se repitan entre corridas (permite comparar, anes usaba 444444)
semillaG2 <- XXXG2 			# Establece el valor de la semilla
fGuardarDatos <- 0			# Decide si guarda (1) o no (0) todos los valores de q, d2 y rho
fGuardarGraficos <- 0			# Decide si guarda (1) o no (0) todos los gráficos en cada ciclo
fGuardarGrafico12Param <- 0		# Decide si guarda (1) o no (0) el grafico de los 12 parámetros
nWrite <- 10000				# Cada cuantos pasos guarda los valores de q, d2 y rho
funcionDens <- "DifCuad"			# [NUEVO] Que funcion usa "DifCuad" = diferencia de cuadrados | "Abs" = valor absoluto


##--- Parámetros a definir para el cálculo teórico mediante diagonalizacion

voltajes <- c((-4:8)*10)   # [mV] c(-20, 0 , 20) Prueba 3 voltajes // c((-4:4)*10) Prueba 9 voltajes // c(5, 25 , 35) Prueba curvas experimentales
T <- 310					# [K] Temperatura (es 310 en el paper del modelo, 298 en experimentos)
p0 <- c(1, 0, 0, 0)			# Estado inicial del sistema
nEstados <- length(p0)			# Addim. Número de estados cinéticos del canal


##--- Parámetros a definir para el algoritmo MCMC

limInf <- c(1e-7, 1e-7, 1e-7, 3.76e-7, 1e-7, 1e-7, 0, 0, 0, 2.07, 0, 0)	# Límite inferior de cada parametro
limSup <- c(1e-1, 1e-1, 1e-1, 3.79e-7, 1e-1, 1e-1, 4, 4, 4, 2.075, 4, 4)	# Límite superior de cada parametro
tipoMesh <- c(0,0,0,0,0,0,1,1,1,1,1,1)	# vector donde se indica si se usa mesh lineal (1) o logaritmico (0)
