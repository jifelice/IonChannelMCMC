# Programa con nucleo de MCMC para ajuste "automatico" de parámetros de modelo
#
# Inputs: q, d2Ini, nPasos, Nini, d2Fin, Kv, limSup, limInf, rho, theorData, expData, sigma
#
# Outputs: q, d, rho

## Algoritmo MCMC

calcularMCMC <- function(q, d2Ini, nPasosBase, Nini, d2Fin, Kv, limSup, limInf, theorData, expData, 
					sigmaIni,told2, nMaxCiclos,factorSigma,factorN,fPasos,cargaMin,cargaMax,
					fGuardarDatos,fGuardarGraficos,fGuardarGrafico12Param,semilla,nWrite,
					ciclos,Nmax,maxCiclosNSigma){ 

	# Esto es por si estamos continuando una corrida (no sigue corriendo)

	if (ciclos == nMaxCiclos){
		print("Se completaron todos los ciclos")
		quit()
		}


	nCoord <- length(q)	# El número de coordenadas sale del largo del vector inicial (q0)

	N <- Nini			# La variable va a ser N, inicialmente toma el valor Nini

	nPasos <- nPasosBase + fPasos*15000*N/10		# Solo para iniciar la matriz MCMC

	MCMC <- matrix(0, ncol=length(q)+2, nrow=nMaxCiclos*nPasosBase)	# Guarda q, d y rho VER EL nrow!!

	MCMCRes <- matrix(0, ncol=length(q)+2, nrow=ceiling(nPasos/nWrite))

	d2despues <- d2Ini

	d2 <- d2Ini				# Para el primer paso	

	av <- c(1:length(q))*0		# vector con el valor de a (para el mesh) para cada parámetro

	sigma <- sigmaIni

	evolucion <- c(1:6)*0				# Guarda valores de ciclo, sigma, d2 y N por ciclo

	indiceN <- 1						# Indice que determina si sigue subiendo N

	indiceSigma <- 1						# Indice que determina que siga bajando N

	fecha <- format(Sys.Date(), "%y%m%d")		# Lee la fecha del sistema y la pasa a "aammdd" (para nombres de archivos)

	fileArranque <- "arranque.txt"			# Nombre del archivo que guarda por si se corta la simulacion

	fileArranqueTT <- "arranqueTT.txt"			# Nombre del archivo que guarda para usar después en el algoritmo de Tratam Termico

	while(d2despues > d2Fin){

		#Hay que calcular nuevamente el av cuando cambia el N (inputs: limInf, limSup, N)
		# Ver si acá o en otro lado para que no lo calcule siempre.
		
		for (imesh in 1:length(q)){
			av[imesh] <- ((limSup[imesh]+Kv[imesh])/(limInf[imesh]+Kv[imesh]))^(1/N)
			}

		nPasos <- nPasosBase + fPasos*15000*N/10		# Ajusta el nro de pasos al N

		iWrite <- 1							# Contador para ir guardando 

		for (i in 1:nPasos){

			qPrevia <- q					# Valor de q previo al cálculo (sería el q_mu)
			d2Previa <- d2					# Valor de d previo al cálculo (sería el d(q_mu))

			coord <- sample(1:nCoord, 1)			# Elige aleatoriamente la coordenada que se mueve
			s <- sample(c(-1,1),1)				# Elige aleatoriamente hacia qué lado se mueve la coordenada
			q[coord] <- (q[coord] + Kv[coord])*(av[coord]^s)-Kv[coord] # Valor de q nuevo (sería el q_nu)
					
			if (q[coord]<limInf[coord]|q[coord]>limSup[coord]|sum(q[7:12])<cargaMin|sum(q[7:12])>cargaMax){		# Prueba que no se vaya de los límites o de la suma de cargas 

				q[coord] <- qPrevia[coord]

				# d2 y rho quedan con el valor que tenían antes

			} else {

				# --------------- Evalúa si se acepta el movimiento ----------------------

				theorData <- calcularptot(voltajes, q, T, p0, tiempos)	# Calcula nuevas corrientes teoricas

				d2Nueva <- dens(theorData, expData)		# Calcula el nuevo valor de d (sería el d(q_nu))

				if (d2Nueva < d2Previa){			# Se acepta porque bajó la distancia
					d2 <- d2Nueva					
				} else {
					deltad2 <- d2Nueva-d2Previa		# es > 1 (porque la d2 no bajó)
					rnd1 <- runif(1)
					pAcept <- exp(-deltad2/(2*sigma^2))
					if (rnd1 <= pAcept) {	# Se acepta
						d2 <- d2Nueva					
					} else {					# Se rechaza 
						q <- qPrevia			# Si se rechaza, q vuelve al valor anterior. No se guarda
						}
					}	# end if d2Nueva < d2Previa


				} #fin if de limites y sumas de cargas


			# ---- Guarda los valores en los vectores (guarda por más que se acepte o no el paso) ---------------------------

			indice <- i + ciclos*nPasos	# Va calculando los pasos 
			MCMC[indice,1:nCoord] <- q
			MCMC[indice,nCoord+1] <- d2

			# ---- Guarda los valores cada cierto tiempo (más adelante los guarda en un .txt) -------------------

			if (i%%nWrite == 0){		
				
				MCMCRes[iWrite,1] <- indice			# Guarda el número de paso
				MCMCRes[iWrite,2:(nCoord+1)] <- q
				MCMCRes[iWrite,nCoord+2] <- d2

				iWrite <- iWrite + 1
				}

			}	# end for de nPasos

			ciclos <- ciclos + 1			# Suma 1 al numero de ciclos despues de hacer nPasos

			# ---- Guarda los valores cada cierto tiempo en un .txt ---------------------------

			nombreArchivoRes <- paste("ValoresRes_ciclo",ciclos,"_",fecha,"_",semilla,".txt",sep="") 
			write.table(MCMCRes, file=nombreArchivoRes, sep = "\t ",col.names=FALSE,row.names=FALSE,quote=FALSE)

			# ---- Guarda los valores al final de un ciclo en un .txt

			pasos <- ciclos*nPasos

			ResultadoCiclo <- t(c(q, d2))

			nombreArchivo <- paste("Valores_ciclo",ciclos,"_",fecha,"_",semilla,".txt",sep="") 
			write.table(ResultadoCiclo, file=nombreArchivo , sep = "\t ",col.names=FALSE,row.names=FALSE,quote=FALSE)


			# ------- Guarda un archivo .jpg con las curvas exp y teor

			if (fGuardarGraficos == 1){

				pasos <- ciclos*nPasos
				nombreArchivo <- paste("curvas_",pasos,"pasos.jpg",sep="") 

				png(filename=nombreArchivo)

				titulo <- paste(pasos, "pasos")

###G  plot(theorData[,1], theorData[,2], main=titulo, col=1, ylim=c(0,1), xlab="Tiempo(ms)", ylab="Po", pch=4)
###G  for (iplot in 3:ncol(theorData)){
###G  points(theorData[,1], theorData[,iplot], col=iplot-1, pch=4)
###G  }
###G  for (iplot in 2:ncol(expData)){
###G  points(expData[,1], expData[,iplot], col=iplot-1)
###G  }
###G  legend(0, 1, legend=c("Experimental", "Teórica"), pch=c(1,4))

				dev.off()
				}


			# ------- Decide como seguir (baja sigma, sube N, etc.) ---------------

			# Calcular d2antes (es la distancia que había antes del nuevo ciclo de nPasos)
			# d2antes es un promedio (salvo en el primer paso que es el valor inicial)

			if (ciclos == 1){
				d2antes <- d2Ini
				}

			# Calcular d2despues (es la distancia que quedo despues del nuevo ciclo de nPasos)
			# d2despues tambien es un promedio

			d2A <- ciclos*nPasos - nPromd2	# Desde el paso A empieza a calcular el promedio de d2despues 
			d2B <- ciclos*nPasos		# Hasta el paso B calcula el promedio de d2despues

			d2despues <- mean(MCMC[d2A:d2B,nCoord+1])			

			evolucion[1] <- ciclos
			evolucion[2] <- d2antes
			evolucion[3] <- d2despues
			evolucion[4] <- (d2antes-d2despues)/d2antes
			evolucion[5] <- sigma
			evolucion[6] <- N

				# --- Guarda la evolución de un ciclo 

				filenameEvolucion <- paste("evolucion_nPasos",nPasos,".txt",sep="_")
				cat(t(evolucion), "\n", file=filenameEvolucion, append=TRUE, sep="\t")


				# ---- Crea el archivo arranque.txt por si se corta la simulacion

				dataCiclo <- paste("ciclos <-", ciclos)
				dataSigma <- paste("sigmaIni <-", sigma)
				datad2 <- paste("d2Ini <-", d2)
				dataN <- paste("Nini <-", N)
				datad2antes <- paste("d2antes <-", d2despues)		#la d2antes del prox paso es la d2antes de este

				cat(dataCiclo, "\n", file=fileArranque, append=FALSE, sep="")
				cat("q <- c(",q[1],",",q[2],",",q[3],",",q[4],",",q[5],",",q[6],",",q[7],",",q[8],",",q[9],",",q[10],",",q[11],",",q[12],")","\n", file=fileArranque, append=TRUE, sep="")
				cat(dataSigma, "\n", file=fileArranque, append=TRUE, sep="")
				cat(datad2, "\n", file=fileArranque, append=TRUE, sep="")
				cat(dataN, "\n", file=fileArranque, append=TRUE, sep="")
				cat(datad2antes, "\n", file=fileArranque, append=TRUE, sep="")


				if ((d2antes-d2despues)/d2antes > told2){		# la distancia(%) bajó más que told2

					sigma <- sigma/factorSigma			# Baja sigma
					indiceN <- 1					# Lo lleva a 1, primer ciclo con ese sigma y N
					indiceSigma <- 1					# Lo lleva a 1, primer ciclo con ese sigma y N

				} else {							# la distancia(%) NO bajó más que told2
				
					if (indiceN == 1 & N <= Nmax & indiceSigma == 1){		# Sube el N una vez

						N <- N*factorN				# Tambien hay que aumentar indiceN, pero lo hace despues

						}

					if (indiceN <= maxCiclosNSigma){		# numero de ciclos con ese sigma y N es menor que el maxCiclosNSigma

						indiceN <- indiceN+1			# sigue con ese sigma y N (solo incrementa indiceN)
                                          }else{
					     sigma <- sigma/factorSigma			# Baja sigma
					     indiceN <- 1
					     indiceSigma <- 1
					}
					}		# end (d2antes-d2despues)/d2antes > told2


			d2antes <- d2despues			# Deja el nuevo valor para comparar con el proximo ciclo

		if (ciclos == nMaxCiclos){
			break
			}

		}	# end while


	# ---

	if (fGuardarDatos == 1){
		filenametxt <- paste("ajuste-corrientes-",fecha,"-",semilla,".txt", sep="")	# Nombre del archivo a guardar

		write.table(MCMC, file=filenametxt, sep = "\t ",col.names=FALSE,row.names=FALSE,quote=FALSE)
		}

	if (fGuardarGrafico12Param == 1){
		#############- GRÁFICO DE PARÁMETROS A LO LARGO DE LOS PASOS (a un .png)

		final <- length(MCMC[,1])

		fecha <- format(Sys.Date(), "%y%m%d")		# Lee la fecha del sistema y la pasa a "aammdd"
 
		paso <- c(1:length(MCMC[1:final,1]))
		parametro <- c("alpha_0_1", "alpha_0_2", "alpha_0_3", "beta_0_1", "beta_0_2",
					"beta_0_3", "zeta_alpha_1", "zeta_alpha_2", "zeta_alpha_3",
					"zeta_beta_1", "zeta_beta_2", "zeta_beta_3")

		q0Villalba <- c(4.5940e-4, 1.3196e-3, 1.5538e-2, 3.7759e-7, 4.4267e-3, 8.6658e-4, 	# Valor inicial de coordenadas (para el HVCN1 son 12 parámetros = 12 coordenadas)
						7.4982e-1, 1.0120, 1.3479e-5, 2.0726, 2.0705, 1.6309e-7)

		filenamepng <- paste("12param-",fecha,"-",semilla,".png", sep="")	# Nombre del archivo a guardar

		png(filename=filenamepng, width = 900, height = 600, unit = "px")

		par(mfrow=c(3,4), cex.axis=1.4, cex.lab=1.75)			# Organizados en 3 filas de 4 columnas cada una

###G for (iParam in 1:6){
		#	plot(paso, MCMC[1:final, iParam], xlab="Número paso", ylab=parametro[iParam], log="y")
###G plot(paso, MCMC[1:final, iParam], xlab="Número paso", ylab=parametro[iParam], ylim=c(limInf[iParam],limSup[iParam]), log="y")
###G abline(h=q0Villalba[iParam], col="red")
###G }

###G for (iParam in 7:12){
		#	plot(paso, MCMC[1:final, iParam], xlab="Número paso", ylab=parametro[iParam])
###G plot(paso, MCMC[1:final, iParam], xlab="Número paso", ylab=parametro[iParam], ylim=c(limInf[iParam],limSup[iParam]))
###G abline(h=q0Villalba[iParam], col="red")
###G }

		dev.off()	

		}

	# --- Guarda el valor de semilla y ÚLTIMOS valores de q, d2 y ciclo 

#	filenametxtfinal <- paste("ResultadosFinal-",semilla,"-",fecha,".txt", sep="")	# Nombre del archivo a guardar
	filenametxtfinal <- paste("ResultadosFinal.txt")				# Nombre del archivo a guardar (sin fecha ni semilla)

	ValoresFinales <- t(c(semilla, q, d2, ciclos))

	write.table(ValoresFinales, file=filenametxtfinal, sep = "\t ",col.names=FALSE,row.names=FALSE,quote=FALSE)


	# ---- Crea el archivo arranqueTT.txt por si se ejecuta después el algoritmo de Tratamiento Termico

	cat("q <- c(",q[1],",",q[2],",",q[3],",",q[4],",",q[5],",",q[6],",",q[7],",",q[8],",",q[9],",",q[10],",",q[11],",",q[12],")","\n", file=fileArranqueTT, append=TRUE, sep="")
	cat("d2Ini <- ", d2, "\n", file=fileArranqueTT, append=TRUE, sep="")
		
	return(MCMC)

	} # end function

