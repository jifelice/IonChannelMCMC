# Guarda los parámetros iniciales de la corrida en un .txt
# Va a guardar todo lo que figura en "nombres"


informeInicialTT <- function (d2Ini,fileExpData,q,cargaMin,cargaMax,nWrite,informeTxtTT,sigmasTT,NsTT,pasosTT){

	datosTT <- matrix(0, ncol=13, nrow=9)
	
	nombresTT <- c("d2Ini","cargaMin","cargaMax","fileExpData","nWrite","q","sigmasTT","NsTT","pasosTT")

	valoresTT <- c(d2Ini,cargaMin,cargaMax,fileExpData,nWrite)

	for (d in 1:length(nombresTT)){
		datosTT[d,1] <- nombresTT[d]
		datosTT[d,2] <- valoresTT[d]
		}

	nCiclosTT <- length(sigmasTT)

###GdatosTT[6,2:13] <- q
###GdatosTT[7,2:(nCiclosTT+1)] <- sigmasTT
###GdatosTT[8,2:(nCiclosTT+1)] <- NsTT
###GdatosTT[9,2:(nCiclosTT+1)] <- pasosTT

	filenametxt <- paste(informeTxtTT,".txt", sep="")
	write.table(datosTT, file=filenametxt, sep = "\t ",col.names=FALSE,row.names=FALSE,quote=FALSE)

	} # end function
