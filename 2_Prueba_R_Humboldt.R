
rm(list=ls())

## Instalar paquetes requeridos
#install.packages("raster")
#install.packages("tera")
#install.packages("spThin")
#install.packages("dismo")
#install.packages("sf")
#install.packages("ENMeval")
#install.packages("rgbif")
#install.packages("dplyr")
#install.packages("readr")
#install.packages("sp")
#install.packages("rgdal")

## Cargar paquetes
library(raster)
library(terra)library(spThin)library(dismo)library(sf)
library(ENMeval)library(rgbif)
library(dplyr)
library(readr)
library("sp")
library("rgdal")

### Cargar capas ambientales. Se requiere un set de 19 capas bioclim con cobertura global. 

#Defnir escritorio donde se encuentran la capas
setwd('~/WorldClim/World')

# Lineas 37-48, se leen las capas y se genera un stack de éstas.
capas=list.files('~/WorldClim/World', pattern='.tif')

lista_capas=list()

for (i in 1:length(capas))

{
	biocapa=raster(capas[i])
	lista_capas[[i]]<-biocapa
}

bio=stack(lista_capas[[1]], lista_capas[[2]], lista_capas[[3]], lista_capas[[4]], lista_capas[[5]], lista_capas[[6]], lista_capas[[7]], lista_capas[[8]], lista_capas[[9]], lista_capas[[10]], lista_capas[[11]], lista_capas[[12]], lista_capas[[13]], lista_capas[[14]], lista_capas[[15]], lista_capas[[16]], lista_capas[[17]], lista_capas[[18]], lista_capas[[19]])

#####

#Lineas 54-65 corresponde a una sección de código tomada del paquete 'rgbif' para descargar una lista de nombre de especies de plantas (https://data-blog.gbif.org/post/downloading-long-species-lists-on-gbif/)

file_url <- "https://data-blog.gbif.org/post/2019-07-11-downloading-long-species-lists-on-gbif_files/global_tree_search_trees_1_3.csv"

long_checklist <- readr::read_tsv(file_url)

# match the names 
gbif_taxon_keys <- 
readr::read_csv(file_url) %>%
head(1000) %>% # use only first 1000 names for testing
pull("Taxon name") %>% # use fewer names if you want to just test 
name_backbone_checklist()  %>% # match to backbone
filter(!matchType == "NONE") %>% # get matched names
pull(usageKey) # get the gbif taxonkeys


# Se genera un pdf vacío donde se guardaran los mapas de distribución para cada especie
pdf('Mapas_Resultados.pdf')

# Loop que permite el proceso iterativo para las 100 especies rqueridas en la prueba.
for (i in 1:100)

{
	# En la primera sección de este loop se descargan los registros para cada especie, filtrando aquellos que solo tiene coordenadas, y eliminando aquellos que están duplicados. Además se eliminan aquellos que caen por fuera de la masa terrestre.
	
	Sp=gbif_taxon_keys[i]
	
	occurrences=occ_data(taxonKey = Sp, hasCoordinate = TRUE,  limit=100, occurrenceStatus = "PRESENT") ##	occurrences=as.data.frame(occurrences$data)	occurrences=occurrences[,c('species', 'decimalLatitude', 'decimalLongitude', 'identifier', 'recordedBy', 'gbifID', 'occurrenceStatus', 'country', 'locality')]
	
	#Remover duplicados
	occurrences <- occurrences[!duplicated(occurrences[, 2:3]), ]
	
	occurrences_xy=occurrences[,c('decimalLongitude', 'decimalLatitude')]
	colnames(occurrences_xy)=c('longitude', 'latitude')
	
	occurrences_xy <- cbind(occurrences_xy, as.data.frame(raster::extract(bio, occurrences_xy)))
	
	##Remover NAs
	occurrences_xy=na.omit(occurrences_xy)
		
	occurrences_xy=data.frame(Sp=unique(occurrences$species), occurrences_xy)
	
	# Se eliminan registros que esten a menos de 5 km
	
	cleaned_occurrences_xy=thin(occurrences_xy, lat.col="latitude", long.col="longitude", spec.col="Sp", reps=1, thin.par=5, locs.thinned.list.return = TRUE, write.files=FALSE, write.log.file=FALSE)
	## Dado que la función "thin" trabaja con un muestreo aleatorio, este proceso debería repetirse varias veces modificando el parametro "reps". Sin embargo, dado que este código está dirigido a especies con pocos registros, para efectos prácticos dejo 1 por defecto.
	
	
	occurrences_xy=occurrences_xy[as.numeric(rownames(cleaned_occurrences_xy[[1]])),]
	
	
	## En este punto, el código define si la especie tiene entre 5 y 15 coordenadas, según lo requerido en la prueba. Si esta condición no se cumple, el proceso se detiene y continua con la siguiente especie.
	
	if(nrow(occurrences_xy)>=5 & nrow(occurrences_xy)<16)
	{
	
	# Coordenadas que definen el "extent" del total de puntos para general los puntos background aleatoreos y generar el modelo de nicho.
	
	maxlat=max(occurrences_xy$latitude)
	minlat=min(occurrences_xy$latitude)
	maxlon=max(occurrences_xy$longitude)
	minlon=min(occurrences_xy$longitude)
	
	coor=occurrences_xy[,2:3]
	coordinates(coor)=c('longitude', 'latitude')
	
	coor=sf::st_as_sf(coor, coords = c('longitude', 'latitude'), crs = raster::crs(bio))
	st_crs(coor) <- '+proj=longlat +datum=WGS84 +no_defs' 
	
	
	### Definición del background
	
	# Para esta prueba, se escogió el método de buffer por punto como área potencial de distribucón de la especies, como también para la selección de puntos aleatoreos. Los buffer fueron definidos a 5 grados por cada punto.
	coor <- sf::st_buffer(coor, dist = 500000) %>% 
	sf::st_union()	%>% 
	sf::st_sf()
	
	
	#Corte de la capas de acuerdo al nuevo background
	bio2=crop(bio2, c(minlon-2, maxlon+2, minlat-2, maxlat+2))
	bio2=mask(bio2, coor)
	
	#Generación de las pseudo-ausencias
	template=as.data.frame(bio2[[1]], xy=TRUE)
	colnames(template)=c('longitude', 'latitude', colnames(template)[3])
	rand=sample(rownames(template), 100000)
	background=template[rand,]
	
	### Partición de los datos para la evaluación de modelos.
	
	### Dado que son especies con pocos datos, se escogió el método Jackknife (ver Shcheglovitova & Anderson 2013)
	
	group.data <- ENMeval::get.jackknife(occurrences_xy[,2:3], background[,1:2])
	
	occs.xy <- occurrences_xy %>% dplyr::select(longitude, latitude)
    
    bg.xy <- background %>% dplyr::select(longitude, latitude)
    
    ### Generación del modelo usando Maxent, y los modelos lineales y cuadráticos, que son los más simples y con menos parámetros dado los pocos puntos. Estos dos modelos se corrieron con tres combinaciones valores de regularización, para un total de seis modelos
    
    Model <- ENMeval::ENMevaluate(occs = occs.xy, 
        bg = bg.xy, partitions = "jackknife", user.grp = group.data, 
        envs = bio2, 
        algorithm = 'maxnet',  tune.args = list(fc = c("L","Q"), rm = 1:3), doClamp=TRUE, parallel=TRUE, numCores=7, parallelType = "doSNOW") #lineal o quadratico, mas simple por los pocos puntos.
        
      results <- eval.results(Model)
        
        ### Selección del modelo basado en el criterio de AIC. Aunque valores de AIC no son recomendados para selección de modelos (Velasco & González-Salazar 2019), por simplicidad de la prueba se usa únicamente este.
        
        bestModel <- results %>% filter(delta.AICc == 0)
        
        pred.seq <- eval.predictions(Model)[[opt.seq$tune.args]]
        plot(bestModel)
        
        Mapa_Modelo_Final=eval.predictions(Model)[[bestModel[1,3]]]
        
        Mapa_Modelo_Final[which(values(Mapa_Modelo_Final)<0.1)]<-NA
        
        #Generación del mapa que será guardado en archivo pdf para todas las especies
        plot(Mapa_Modelo_Final, main=unique(occurrences$species))
        plot(coor, add=T)
             
	}
	
	
	if(nrow(occurrences_xy)<5 | nrow(occurrences_xy)>15) {next}
	
	
}

dev.off()



