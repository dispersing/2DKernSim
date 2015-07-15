##################################
#     Two-dimentional kernel     #
#  resampling of dispersal data  #
#                                #
#             INDEX              #
#      -1. load library(ies)     #
#       0. assign a dataset      #
#       1. randomize theta,      #
#  create 2-dimentional kernel,  #
#          and resample          #
#    2. find distribution mean   #
#       and standard error       #
#         3. plot data           #
#                                #
#       ~written by isadore nabi~#
##################################

# -1. load library(ies)
     library(MASS) #for kde2d and plotting functions

# 0. assign a dataset
	# 0.0 clear objects and set seed
	rm(list = ls())
	set.seed(6174) #Kaprekar's constant

	# 0.1 100 data points with x and y coordinates
	data.points <- 100
	bi.var.norm <- cbind(rnorm(n = data.points, mean = 0 , sd = 15 ) , rnorm(n = data.points, mean = 0 , sd = 15 ))
	disp.data <- (bi.var.norm[,1]^2 + bi.var.norm[,2]^2)^0.5

	# 0.2 plot the data and histogram (distance distribution)
	par(mfrow=c(2,2) , pin = c(2,2))
	plot(bi.var.norm , type = "n")
	abline(h = 0 , v = 0 , lwd = 2 , col = "grey70" , lty = 2)
	points(bi.var.norm , pch = 16 , col = "blue")
	arrows(rep(0 , data.points) , rep(0 , data.points) , bi.var.norm[,1] , bi.var.norm[,2] , length = 0 , col = "red")
	hist(disp.data , main = "")

# 1. randomize theta, create 2-dimentional kernel, and resample
	# 1.0 loop
	samples <- 250 #remember that it is going to be 4 per sample because we're extracing the Cartesian coordinates
	breaks <- max(ceiling(disp.data)) # finds largest distance and creates that many breaks to correspond to one per unit distance

	w.mat <- matrix(data = NA , nrow = samples*4 , ncol = breaks , byrow = T) # for density distribution data (as omega) pre-allocation
	f.mat <- matrix(data = NA , nrow = samples*4 , ncol = breaks , byrow = T) # for distance distribution data (as f) pre-allocation

	distances <- 1:ncol(w.mat)

	for (i in 0:(samples-1)){
		# 1.1 draw theta for number of dispersal sampels
		theta <- round(runif(length(disp.data) , max = 360) , 1)
	
		# 1.2 convert polar cordinates to bearing
		polar2cart <- function(x , y , dist , bearing , as.deg = FALSE){
			  if(as.deg){
			    bearing = bearing * pi / 180
			  }
		  new.x <- x + dist*sin(bearing)
		  new.y <- y + dist*cos(bearing)
		  return(list("x" = new.x,"y" = new.y))
		}
		points.2D <- polar2cart(0 , 0 , disp.data, theta)
	
		# 1.3 create 2D pdf
		kernel.2D <- kde2d(points.2D$x , points.2D$y , n = breaks*2)

		# 1.4 z-values for a sample from pole across a Cartetian ray
		max.x <- sort(kernel.2D$z[(breaks+1):(breaks*2) , breaks] , decreasing = T)
		min.x <- sort(kernel.2D$z[1:breaks , breaks] , decreasing = T)
		max.y <- sort(kernel.2D$z[breaks , 1:breaks] , decreasing = T)
		min.y <- sort(kernel.2D$z[breaks , (breaks+1):(breaks*2)] , decreasing = T)

		w.mat[(i*4+1),] <- max.x/sum(max.x)
		w.mat[(i*4+2),] <- min.x/sum(min.x)
		w.mat[(i*4+3),] <- max.y/sum(max.y)
		w.mat[(i*4+4),] <- min.y/sum(min.y)

		f.r.max.x <- max.x*2*pi*distances
		f.r.min.x <- min.x*2*pi*distances
		f.r.max.y <- max.y*2*pi*distances
		f.r.min.y <- min.y*2*pi*distances

		f.mat[(i*4+1),] <- f.r.max.x/sum(f.r.max.x)
		f.mat[(i*4+2),] <- f.r.min.x/sum(f.r.min.x)
		f.mat[(i*4+3),] <- f.r.max.y/sum(f.r.max.y)
		f.mat[(i*4+4),] <- f.r.min.y/sum(f.r.min.y)
	}

	# 1.4 structural checks
		#1.4.0 head, tails, and row sums
		w.mat[1:5 , 1:5] # head should be larger than the tail (next line below)
		w.mat[1:5 , (ncol(w.mat)-5):ncol(w.mat)] # tail
		w.mat[(nrow(w.mat)-5):nrow(w.mat) , 1:5] # head should be larger than the tail (next line below)
		w.mat[(nrow(w.mat)-5):nrow(w.mat) , (ncol(w.mat)-5):ncol(w.mat)] # tail
		apply(w.mat , 1 , sum) # summing to unity means the sum of rows all should = 1
	
		f.mat[1:5 , 1:5] # head should be larger than the tail (next line below)
		f.mat[1:5 , (ncol(f.mat)-5):ncol(f.mat)] # tail
		f.mat[(nrow(f.mat)-5):nrow(f.mat) , 1:5] # head should be larger than the tail (next line below)
		f.mat[(nrow(f.mat)-5):nrow(f.mat) , (ncol(f.mat)-5):ncol(f.mat)] # tail
		apply(f.mat , 1 , sum) # summing to unity means the sum of rows all should = 1

		#1.4.1 three example plots
		par(mfrow = c(1,3) , pin = c( 2,2))
		plot(points.2D , pch = 21 , bg = "#50505050") # points quadri-sected
		abline(h = 0 , v = 0 , col = "red" , lty = 2)
		contour(kernel.2D) # contour plot with example segment
		arrows(0,0,max(kernel.2D$x),0 , length = 0 , col = "red" , lwd = 2 , lty = 2)
		persp(kernel.2D , phi = 15 , theta = 235 , xlab = "Lateral axis" , ylab = "Longitudinal axis" , zlab = "Density" , shade = 0.9 , expand = 0.9 , col = "white" , lwd = 0.5) # 3D plot, because you have to

# 2. find distribution mean and standard error
	# 2.0 density distance pdf
	wr.mean <- apply(w.mat , 2 , mean)
	wr.sem <- apply(w.mat , 2 , sd)/(nrow(w.mat)^0.5)

	# 2.1 location distance pdf
	fr.mean <- apply(f.mat , 2 , mean)
	fr.sem <- apply(f.mat , 2 , sd)/(nrow(f.mat)^0.5)

# 3. plot data
	# Discrete box plots
	par(mfrow = c(1,2) , mar = c(3,3,0,0) , oma = c(1,1,1,1))
	# density distribution
	boxplot(w.mat , medcol = NA , axes = F , outpch = 3 , outcol = "blue" , outcex = 0.5 , ylim = c(0, max(w.mat)) , col = "grey95" , whisklty = 1)
	abline(h = 0 , lwd = 0.75 , col = "black")
	segments(distances, wr.mean - wr.sem, distances, wr.mean + wr.sem , lwd = 3 , col = "red")
	box()
	axis(1)
	axis(2 , las = 1)
	mtext("Distance from pole" , side = 1 , line = 2.25)
	mtext("Probability" , side = 2 , line = 3)

	# distance distribution
	boxplot(f.mat , medcol = NA , axes = F , outpch = 3 , outcol = "blue" , outcex = 0.5 , ylim = c(0,max(f.mat)) , col = "grey95" , whisklty = 1)
	abline(h = 0 , lwd = 0.75 , col = "black")
	segments(distances, fr.mean - fr.sem, distances, fr.mean + fr.sem , lwd = 3 , col = "red")
	box()
	axis(1)
	axis(2 , las = 1)
	mtext("Distance from pole" , side = 1 , line = 2.25)

	# Continuous box plots
	box.stats.w <- boxplot(w.mat , plot = F)$stats[-3,]
	box.stats.f <- boxplot(f.mat , plot = F)$stats[-3,]

	par(mfrow = c(1,2) , mar = c(3,3,0,0) , oma = c(1,1,1,1))
	plot(0 , type = "n" , ylim = c(0, max(w.mat)) , xlim = c(0,length(box.stats.w[4,])) , xlab = "Distance from pole" , ylab = "" , las = 1)
	polygon(c(1:length(box.stats.w[4,]),length(box.stats.w[4,]):1),c(box.stats.w[4,],rev(box.stats.w[1,])) , col = rgb(0,0,1,1/3) , border = F)
	polygon(c(1:length(box.stats.w[4,]),length(box.stats.w[4,]):1),c(box.stats.w[3,],rev(box.stats.w[2,])) , col = rgb(0,0,.75,2/3) , border = F)
	lines(wr.mean , col = "red" , lwd = 2 , lend = 1)
	abline(h = 0 , lwd = 0.75 , col = "black")
	box()
	axis(1 , labels = F)
	mtext("Probability" , side = 2 , line = 3)
	mtext("Distance from pole" , side = 1 , line = 2.25)

	plot(0 , type = "n" , ylim = c(0, max(f.mat)) , xlim = c(0, length(box.stats.f[4,])) , xlab = "Distance from pole" , ylab = "" , las = 1)
	polygon(c(1:length(box.stats.f[4,]),length(box.stats.f[4,]):1),c(box.stats.f[4,],rev(box.stats.f[1,])) , col = rgb(0,0,1,1/3) , border = F)
	polygon(c(1:length(box.stats.f[4,]),length(box.stats.f[4,]):1),c(box.stats.f[3,],rev(box.stats.f[2,])) , col = rgb(0,0,.75,2/3) , border = F)
	lines(fr.mean , col = "red" , lwd = 2 , lend = 1)
	abline(h = 0 , lwd = 0.75 , col = "black")
	box()
	axis(1 , labels = F)
	mtext("Distance from pole" , side = 1 , line = 2.25)
