# Guarda los parámetros iniciales de la corrida en un .txt
# Va a guardar todo lo que figura en "nombres"


informeInicial <- function (d2Ini,d2Fin,nPasos,nMaxCiclos,told2,nPromd2,Nini,fileExpData,sigmaIni,semilla,q,
					factorSigma,factorN,cargaMin,cargaMax,nWrite,informeTxt,Nmax,maxCiclosNSigma){

	datos <- matrix(0, ncol=13, nrow=18)
	
	nombres <- c("d2Ini","d2Fin","nPasos","nMaxCiclos","cargaMin","cargaMax","told2","nPromd2","Nini","Nmax",
				"maxCiclosNSigma","fileExpData","sigmaIni", "semilla","factorSigma","factorN","nWrite","q")


	valores <- c(d2Ini,d2Fin,nPasos,nMaxCiclos,cargaMin,cargaMax,told2,nPromd2,Nini,Nmax,maxCiclosNSigma,
				fileExpData,sigmaIni,semilla,factorSigma,factorN,nWrite,q)


	for (d in 1:length(nombres)){
		datos[d,1] <- nombres[d]
		datos[d,2] <- valores[d]
		}

	datos[18,2:13] <- q

	
	filenametxt <- paste("informeInicial-", semilla, ".txt", sep="")
	write.table(datos, file=filenametxt, sep = "\t ",col.names=FALSE,row.names=FALSE,quote=FALSE)

	} # end function
