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
#directorioFuentes <- "D:/20220414/1/"	# Donde estan los archivos fuente

#############- Leer parametros:

source("Parametros.R", echo=TRUE)							# Archivo con parametros iniciales (antes del primer codigo)
#source(paste(directorioFuentes, "Parametros.R", sep=""), echo=TRUE)	# Archivo con parametros iniciales (antes del primer codigo)
source("../ParametrosTT.R", echo=TRUE)							# Archivo con valores de sigma, N y numero de pasos
#source(paste(directorioFuentes, "ParametrosTT.R", sep=""), echo=TRUE)	# Archivo con valores de sigma, N y numero de pasos


#############- Cargo las funciones:

source(paste(directorioFuentes, "AJUSTEcoreDIAG.R", sep=""))	# Establece en que archivo esta el nucleo para diagonalizar
source(paste(directorioFuentes, "AJUSTEmatrizHVCN1.R", sep=""))	# Establece en que archivo esta el calculo de la matriz W

if (funcionDens == "DifCuad"){						# [NUEVO] Definir la densidad de probabilidad 
	source(paste(directorioFuentes, "AJUSTEdensfunct.R", sep=""))
	source(paste(directorioFuentes, "AJUSTEcoreMCMC-TT.R", sep=""))	# Establece en que archivo esta el nucleo el MCMC
	} else if (funcionDens == "Abs") {
		source(paste(directorioFuentes, "AJUSTEdensfunctAbs.R", sep=""))		
		source(paste(directorioFuentes, "AJUSTEcoreMCMCAbs.R", sep=""))	# Establece en que archivo esta el nucleo el MCMC
		} else {
		print("la funcion ingresada no es correcta")
		}

source(paste(directorioFuentes, "AJUSTEinformeTT.R", sep=""))		# En este archivo se arma un .txt con los valores iniciales


#############- Leer valores experimentales de corrientes:

archivoExp <- "expData.txt"			# Nombre del archivo con los datos experimentales

fileExpData <- paste("../", archivoExp, sep="")	# otros: "expDataPrueba3V.txt", "expDataPrueba9V.txt", "210304_002.IV IO.1(5,25,35mV)Res.txt"
expData <- read.table(fileExpData, header=FALSE)
tiempos <- expData[,1]				# Tiempos a calcular para las curvas teóricas


#############- Para el algoritmo MCMC:

Kv <- c(1:length(q))*0		# Inicia el vector Kv con el valor de K para cada parámetro

for (iTipo in 1:length(tipoMesh)){
	Kv[iTipo] <- limSup[iTipo]*100*tipoMesh[iTipo]	# K = 0 (escala logarítmica) y K >> q max (escala lineal)
	}

##############- Inicio desde el lugar donde finalizo el algoritmo anterior o donde se corto el de Tratamiento Termico

primerCiclo <- 1

archivos <- list.files()		# guarda en ese vector el nombre de los archivos del directorio

	if ("arranqueTTreanuda.txt" %in% archivos == TRUE){ # busca el archivo en la lista

		source("arranqueTTreanuda.txt", echo=TRUE)		# Carga el archivo con q y d2 para continuar el Tratamiento Termico

		informeTxtTT <- "informeReanudacionTT"

		theorData <- calcularptot(voltajes, q, T, p0, tiempos)	# Calcula con el q que termino

	} else {

		if ("arranqueTT.txt" %in% archivos == TRUE){	# busca el archivo en la lista
	
			source("arranqueTT.txt", echo=TRUE)		# Carga el archivo con q y d2 para empezar el Tratamiento Termico

			informeTxtTT <- "informeInicioTT"

			theorData <- calcularptot(voltajes, q, T, p0, tiempos)	# Calcula con el q que termino
	
		} else {		# Si no está va a inicio aleatorio

			print("No hay un archivo de arranqueTT")
			q()									# Si no hay archivo de arranque, sale del programa

		}

	}


#############- Guardar parámetros iniciales/al reanudar en archivo de texto

informeInicialTT(d2Ini,fileExpData,q,cargaMin,cargaMax,nWrite,informeTxtTT,sigmasTT,NsTT,pasosTT)

#############- Ajuste de parámetros de q mediante MCMC

ajuste <- calcularMCMCTT(Kv, limSup, limInf, sigmasTT, NsTT, pasosTT, cargaMin, cargaMax, expData,q,d2Ini,
					fGuardarDatos,fGuardarGraficos,primerCiclo)


end_time <- Sys.time()					# Finaliza el timer
print(end_time-start_time)


#######################################################################################################
#																	#
#								FIN									#
#																	#
#######################################################################################################
