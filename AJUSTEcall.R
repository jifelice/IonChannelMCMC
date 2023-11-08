# Lee los datos de una corriente experimental y por otro lado va ajustando una
# corriente teórica, para ajustar los parámetros y llegar a los que hagan coincidir 
# las dos curvas. Lo hace en forma automática hasta que se cumple una cierta distancia
# entre las curvas (iniciado en agosto de 2021)

# 1. Leer los valores de corrientes experimentales.
# 2. Definir condiciones iniciales del vector q y voltajes (V)
# 2. Para ese q y cada valor de V, calculará una matriz W y calculará la Po (p_open) teórica.
# 3. Calcula la función densidad de probabilidad (comparando valores experimentales y valores teóricos).
# 4. Aleatoriamente cambia una de las coordenadas de q (es uno de los parámetros y vuelve a calcular la función densidad
#	de probabilidad. 

########################################################################################################################


start_time <- Sys.time()					# Inicia el timer

directorioFuentes <- "/home/corridas/MCMC-VG/fuentes/"	# Donde estan los archivos fuente

#############- Leer parametros:

source("Parametros.R", echo=TRUE)							# Archivo con parametros iniciales
#source(paste(directorioFuentes, "Parametros.R", sep=""), echo=TRUE)	# Archivo con parametros iniciales

#############- Cargo las funciones:

source(paste(directorioFuentes, "AJUSTEcoreDIAG.R", sep=""))	# Establece en que archivo esta el nucleo para diagonalizar
source(paste(directorioFuentes, "AJUSTEmatrizHVCN1.R", sep=""))	# Establece en que archivo esta el calculo de la matriz W

if (funcionDens == "DifCuad"){						# [NUEVO] Definir la densidad de probabilidad 
	source(paste(directorioFuentes, "AJUSTEdensfunct.R", sep=""))
	source(paste(directorioFuentes, "AJUSTEcoreMCMC.R", sep=""))	# Establece en que archivo esta el nucleo el MCMC
	} else if (funcionDens == "Abs") {
		source(paste(directorioFuentes, "AJUSTEdensfunctAbs.R", sep=""))		
		source(paste(directorioFuentes, "AJUSTEcoreMCMCAbs.R", sep=""))	# Establece en que archivo esta el nucleo el MCMC
		} else {
		print("la funcion ingresada no es correcta")
		}

source(paste(directorioFuentes, "AJUSTEinforme.R", sep=""))		# En este archivo se arma un .txt con los valores iniciales


#############- Leer valores experimentales de corrientes:

archivoExp <- "expData.txt"			# Nombre del archivo con los datos experimentales

fileExpData <- paste("../", archivoExp, sep="")	# otros: "expDataPrueba3V.txt", "expDataPrueba9V.txt", "210304_002.IV IO.1(5,25,35mV)Res.txt"
expData <- read.table(fileExpData, header=FALSE)
tiempos <- expData[,1]				# Tiempos a calcular para las curvas teóricas


#############- Para el algoritmo MCMC:

q <- c(1:12)*0			# Inicia el vector q
Kv <- c(1:length(q))*0		# Inicia el vector Kv con el valor de K para cada parámetro

for (iTipo in 1:length(tipoMesh)){
	Kv[iTipo] <- limSup[iTipo]*100*tipoMesh[iTipo]	# K = 0 (escala logarítmica) y K >> q max (escala lineal)
	}

##############- Inicio aleatorio o inicio desde el lugar donde se corto

archivos <- list.files()		# guarda en ese vector el nombre de los archivos del directorio

if ("arranque.txt" %in% archivos == TRUE){	# busca el archivo en la lista

	source("arranque.txt", echo=TRUE)		# Carga el archivo con parametros para empezar

	informeTxt <- "informeReanudacion-"

	theorData <- calcularptot(voltajes, q, T, p0, tiempos)	# Calcula con el q que termino
	
	} else {		# Si no está va a inicio aleatorio

	ciclos <- 0		# para ir contando el número de ciclos de nPasos que va haciendo

	factorInicio <- 0.95			# para evitar los bordes

	while (sum(q[7:12])<= cargaMin | sum(q[7:12])>cargaMax){

		rndIniVect <- runif(12)			# elige un aleatorio entre 0 y 1
		for (iIni in 1:length(limInf)){
			q[iIni] <- limInf[iIni]+factorInicio*rndIniVect[iIni]*(limSup[iIni]-limInf[iIni])

			}
		informeTxt <- "informeInicial-"

		}

	theorData <- calcularptot(voltajes, q, T, p0, tiempos)	# Iniciales

	d2Ini <- dens(expData, theorData)		# d inicial con los valores iniciales

	if (funcionDens == "DifCuad") {

		sigmaIni <- sqrt(d2Ini)/3			# Mínimo sigma para que camine independientemente del inicio	
		rhoIni <- exp(-d2Ini/(2*sigmaIni^2))		# rho  inicial con los valores iniciales

	} else if (funcionDens == "Abs") {

		sigmaIni <- d2Ini/3   		# Mínimo sigma para que camine independientemente del inicio	
		rhoIni <- exp(-d2Ini/sigmaIni)							# rho  inicial con los valores iniciales

		}

	}

#############- Guardar parámetros iniciales/al reanudar en archivo de texto

informeInicial(d2Ini,d2Fin,nPasosBase,nMaxCiclos,told2,nPromd2,Nini,fileExpData,sigmaIni,semilla,q,factorSigma,
			factorN,cargaMin,cargaMax,nWrite,informeTxt,Nmax,maxCiclosNSigma)

#############- Ajuste de parámetros de q mediante MCMC

ajuste <- calcularMCMC(q, d2Ini, nPasosBase, Nini, d2Fin, Kv, limSup, limInf, theorData, expData, 
				sigmaIni, told2, nMaxCiclos, factorSigma, factorN,fPasos,cargaMin,cargaMax,
				fGuardarDatos,fGuardarGraficos,fGuardarGrafico12Param,semilla,nWrite,ciclos,
				Nmax,maxCiclosNSigma)


end_time <- Sys.time()					# Finaliza el timer
print(end_time-start_time)



#######################################################################################################
#																	#
#								FIN									#
#																	#
#######################################################################################################
