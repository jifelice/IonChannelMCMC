# Programa con nucleo de MCMC para ajuste mediante "tratamiento termico" de los parámetros de modelo
# En principio, sería para usar después de correr el algoritmo de ajuste automático
#
# Va leer resultados del archivo arranque.Txt y parametros (sigma, N y número de pasos) de ParametrosTT.R


## Algoritmo MCMC

calcularMCMCTT <- function(Kv, limSup, limInf, sigmasTT, NsTT, pasosTT, cargaMin, cargaMax, expData,q,d2Ini,
					fGuardarDatos,fGuardarGraficos, primerCiclo){ 


	nCoordTT <- length(q)	# El número de coordenadas sale del largo del vector inicial (q0)

	nCiclosTT <- length(sigmasTT)

	MCMCTT <- matrix(0, ncol=length(q)+2, nrow=sum(pasosTT))		# Guarda q, d y rho

	MCMCResTT <- matrix(0, ncol=length(q)+2, nrow=ceiling(sum(pasosTT)/nWrite))

	fileArranqueTTreanuda <- "arranqueTTreanuda.txt"	# Nombre del archivo que guarda por si se corta la simulacion

	qPreviaTT <- q						# El valor que entra como input a este archivo

	d2 <- d2Ini							# Para el primer paso	

	d2antes <- d2Ini						# Para el primer paso

	iWrite <- 1

	for (iCicloTT in primerCiclo:nCiclosTT){

		av <- c(1:length(q))*0				# vector con el valor de a (para el mesh) para cada parámetro

		evolucionTT <- c(1:6)*0				# Guarda valores de ciclo, sigma, d2 y N por ciclo

		sigma <- sigmasTT[iCicloTT]			# Indica el valor de sigma para cada ciclo

		N <- NsTT[iCicloTT]				# Indica el valor de N para cada ciclo

		pasos <- pasosTT[iCicloTT]			# Indica el valor de sigma para cada ciclo

		#Hay que calcular nuevamente el av cuando cambia el N (inputs: limInf, limSup, N)
		
		for (imesh in 1:length(q)){
			av[imesh] <- ((limSup[imesh]+Kv[imesh])/(limInf[imesh]+Kv[imesh]))^(1/N)
			}

		for (i in 1:pasos){

			qPreviaTT <- q					# Valor de q previo al cálculo (sería el q_mu)
			d2PreviaTT <- d2					# Valor de d previo al cálculo (sería el d(q_mu))

			# Elige que parametro va a cambiar
			coord <- sample(1:nCoordTT, 1)					# Elige aleatoriamente la coordenada que se mueve
			s <- sample(c(-1,1),1)							# Elige aleatoriamente hacia qué lado se mueve la coordenada
                        rrr=runif(1)
                        if(coord<=6){
			q[coord] <- (q[coord] + Kv[coord])*(av[coord]^(rrr*s))-Kv[coord]	# Valor de q nuevo (sería el q_nu)
                        } else {
			chuqui <- (q[coord] + Kv[coord])*(av[coord]^s)-Kv[coord]
                        dddqqq=chuqui-q[coord]
                        q[coord]=q[coord] + dddqqq*rrr
                        } 
					
			# Se fija si se va de los límites de los parámetros o de la suma de cargas 

			if (q[coord]<limInf[coord]|q[coord]>limSup[coord]|sum(q[7:12])<cargaMin|sum(q[7:12])>cargaMax){		

				q[coord] <- qPreviaTT[coord]		# Vuelve al q anterior (d2 y rho quedan como estaban)

				} else {

				# ----- Para los que no se van de los limites, evalua si se acepta el movimiento ---------------

				theorData <- calcularptot(voltajes, q, T, p0, tiempos)	# Calcula nuevas corrientes teoricas

				d2NuevaTT <- dens(theorData, expData)		# Calcula el nuevo valor de d (sería el d(q_nu))

				if (d2NuevaTT < d2PreviaTT){			# Se acepta porque bajó la distancia

					d2 <- d2NuevaTT					

				} else {

					deltad2 <- d2NuevaTT-d2PreviaTT		# es > 1 (porque la d2 no bajó)

					rnd1 <- runif(1)

					pAcept <- exp(-deltad2/(2*sigma^2))

					if (rnd1 <= pAcept) {			# Se acepta

						d2 <- d2NuevaTT					

					} else {					# Se rechaza 

						q <- qPreviaTT			# Si se rechaza, q vuelve al valor anterior. No se guarda

						}

					}	# end if d2NuevaTT < d2PreviaTT

				} #fin if de limites y sumas de cargas

			# ---- Guarda los valores en los vectores (guarda por más que se acepte o no el paso) -----------------

			indiceTT <- i + sum(pasosTT[1:iCicloTT]) - pasosTT[1]			# Va calculando los pasos 
			MCMCTT[indiceTT,1:nCoordTT] <- q
			MCMCTT[indiceTT,nCoordTT+1] <- d2

			# ---- Guarda los valores cada cierto tiempo (más adelante los guarda en un .txt) -------------------

			if (i%%nWrite == 0){						# Si i es múltiplo de nWrite, guarda		
				
				MCMCResTT[iWrite,1] <- indiceTT			# Guarda el número de paso
				MCMCResTT[iWrite,2:(nCoordTT+1)] <- q
				MCMCResTT[iWrite,nCoordTT+2] <- d2

				iWrite <- iWrite + 1
				}

			}	# end for pasos


			# ---- Guarda los valores cada cierto tiempo en un .txt ---------------------------

			inicioGuardar <- 1+ceiling(sum(pasosTT[1:iCicloTT])/nWrite)-ceiling(sum(pasosTT[1])/nWrite)

			nombreArchivoRes <- paste("ValoresRes_ciclo_TT_",iCicloTT,".txt",sep="") 
			write.table(MCMCResTT[inicioGuardar:(iWrite-1),], file=nombreArchivoRes, sep = "\t ",col.names=FALSE,row.names=FALSE,quote=FALSE)

			# ---- Guarda los valores al final de un ciclo en un .txt

			ResultadoCiclo <- t(c(q, d2))

			nombreArchivo <- paste("Valores_ciclo_TT_",iCicloTT,".txt",sep="") 
			write.table(ResultadoCiclo, file=nombreArchivo , sep = "\t ",col.names=FALSE,row.names=FALSE,quote=FALSE)

			# ------- Guarda un archivo .jpg con las curvas exp y teor

			if (fGuardarGraficos == 1){

				pasosJPG <- sum(pasosTT[1:iCicloTT])
				nombreArchivo <- paste("curvas_",pasosJPG,"pasos.jpg",sep="") 

				png(filename=nombreArchivo)

				titulo <- paste(pasosJPG, "pasos")

###G plot(theorData[,1], theorData[,2], main=titulo, col=1, ylim=c(0,1), xlab="Tiempo(ms)", ylab="Po", pch=4)
###G for (iplot in 3:ncol(theorData)){
###G points(theorData[,1], theorData[,iplot], col=iplot-1, pch=4)
###G }
###G for (iplot in 2:ncol(expData)){
###G points(expData[,1], expData[,iplot], col=iplot-1)
###G }

###G legend(0, 1, legend=c("Experimental", "Teórica"), pch=c(1,4))

				dev.off()
				}

			# Crea el archivo de evolucion 

			pasosEvol <- sum(pasosTT[1:iCicloTT])

			d2A <- pasosEvol-nPromd2		# Desde el paso A empieza a calcular el promedio de d2despues 
			d2B <- pasosEvol 				# Hasta el paso B calcula el promedio de d2despues
	
			d2despues <- mean(MCMCTT[d2A:d2B,nCoordTT+1])			

			evolucionTT[1] <- iCicloTT
			evolucionTT[2] <- d2antes
			evolucionTT[3] <- d2despues
			evolucionTT[4] <- (d2antes-d2despues)/d2antes
			evolucionTT[5] <- sigma
			evolucionTT[6] <- N

			filenameEvolucion <- "evolucionTT.txt"
			cat(t(evolucionTT), "\n", file=filenameEvolucion, append=TRUE, sep="\t")

			# ---- Crea el archivo arranqueTTreanuda.txt por si se ejecuta después el algoritmo de Tratamiento Termico

			cat("primerCiclo <- ", iCicloTT+1, "\n", file=fileArranqueTTreanuda, append=FALSE, sep="")
			cat("q <- c(",q[1],",",q[2],",",q[3],",",q[4],",",q[5],",",q[6],",",q[7],",",q[8],",",q[9],",",q[10],",",q[11],",",q[12],")","\n", file=fileArranqueTTreanuda, append=TRUE, sep="")
			cat("d2Ini <- ", d2, "\n", file=fileArranqueTTreanuda, append=TRUE, sep="")
		
		d2antes <- d2despues			# es solo para el archivo evolucion

		} # end for iCicloTT

	# --- Guarda todos los datos calculados

	if (fGuardarDatos == 1){
		filenametxt <- "ajuste-corrientes.txt"		# Nombre del archivo a guardar

		write.table(MCMCTT, file=filenametxt, sep = "\t ",col.names=FALSE,row.names=FALSE,quote=FALSE)
		}

	if (fGuardarGrafico12Param == 1){
		#############- GRÁFICO DE PARÁMETROS A LO LARGO DE LOS PASOS (a un .png)

		final <- length(MCMCTT[,1])

		fecha <- format(Sys.Date(), "%y%m%d")		# Lee la fecha del sistema y la pasa a "aammdd"
 
		paso <- c(1:length(MCMCTT[1:final,1]))
		parametro <- c("alpha_0_1", "alpha_0_2", "alpha_0_3", "beta_0_1", "beta_0_2",
					"beta_0_3", "zeta_alpha_1", "zeta_alpha_2", "zeta_alpha_3",
					"zeta_beta_1", "zeta_beta_2", "zeta_beta_3")

		filenamepng <- paste("12param-",fecha,".png", sep="")	# Nombre del archivo a guardar

		png(filename=filenamepng, width = 900, height = 600, unit = "px")

		par(mfrow=c(3,4), cex.axis=1.4, cex.lab=1.75)			# Organizados en 3 filas de 4 columnas cada una

###G for (iParam in 1:6){
		#	plot(paso, MCMCTT[1:final, iParam], xlab="Número paso", ylab=parametro[iParam], log="y")
###G plot(paso, MCMCTT[1:final, iParam], xlab="Número paso", ylab=parametro[iParam], ylim=c(limInf[iParam],limSup[iParam]), log="y")
###G }

###G for (iParam in 7:12){
		#	plot(paso, MCMCTT[1:final, iParam], xlab="Número paso", ylab=parametro[iParam])
###G plot(paso, MCMCTT[1:final, iParam], xlab="Número paso", ylab=parametro[iParam], ylim=c(limInf[iParam],limSup[iParam]))
###G }

		dev.off()	

		}

	# --- Guarda el valor de semilla y ÚLTIMOS valores de q, d2 y ciclo 

	filenametxtfinal <- paste("ResultadosFinalTT.txt")				# Nombre del archivo a guardar

	ValoresFinalesTT <- t(c(q, d2))

	write.table(ValoresFinalesTT, file=filenametxtfinal, sep = "\t ",col.names=FALSE,row.names=FALSE,quote=FALSE)
		
	return(MCMCTT)

	} # end function
