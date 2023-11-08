# Definir la función de densidad de probabilidad para ajustar los parámetros
# del modelo. La función es d(q)=sum[Ci,exp-Ci,teor(q)]^2
#
# Inputs: expData (datos t vs I experimentales) y theorData (un valor de I o Po para cada tiempo).


dens <- function (expData, theorData)  {

	d2 <- 0			# Para que empiece en cero

        pesod2=0
	for (itime in 1:length(tiempos)){
		for (icol in 2:ncol(expData)){			# Hay una columna por cada voltaje	
			exp <- expData[itime,icol]
			teor <- theorData[itime,icol]
			d2 <- d2 + (exp-teor)^2
                        pesod2=pesod2+1
			}
		}
	d2 <- d2/pesod2
	return(d2)
	}
