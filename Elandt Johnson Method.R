
#************** Elandt–Johnson Method for converting five-age mortality rates to single-age **********

#*** import lx from life table ***
lx.five<-read.table("D:\\Book1.txt")
lx.five<-as.matrix(lx.five)

lx.single <- rep(999,97)
lx.single[1] <- 100000

coeffs.ages1_9 <- matrix(
c(
1,		0,		0,		0,		0,		0,
0.56203,	0.7176,	-0.4784,	0.283886,	-0.100716,	0.0156,
0.273392,	1.047199,	-0.531911,	0.2992,	-0.103747,	0.015867,
0.096491,	1.1088,	-0.328533,	0.1728,	-0.058358,	0.0088,
0,		1,		0,		0,		0,		0,
-0.041667,	0.798,	0.354667,	-0.152,	0.048,	-0.007,
-0.048872,	0.5616,	0.6656,	-0.240686,	0.072758,	-0.0104,
-0.037281,	0.3332,	0.888533,	-0.2448,	0.070147,	-0.0098,
-0.018379,	0.1408,	1.001244,	-0.160914,	0.043116,	-0.005867),9,6,byrow=T)


coeffs.ages10_74 <- matrix(
c(
0.008064,	-0.07392,	0.88704,	0.22176,	-0.04928,	0.006336,
0.011648,	-0.09984,	0.69888,	0.46592,	-0.08736,	0.010752,
0.010752,	-0.08736,	0.46592,	0.69888,	-0.09984,	0.011648,
0.006336,	-0.04928,	0.22176,	0.88704,	-0.07392,	0.008064),4,6,byrow=T)


lx.single[2:10] <- coeffs.ages1_9 %*% lx.five[2:7]; ## matrix multiplication;


lx.calc <- lx.five[-2];  ## matrix without age group 1-4;
lx.calc <- cbind(lx.calc, seq(0,85,by=5))
count <- 0

for (i in seq(10,70,by=5)){

	lx.single[i+1] <- lx.calc[which(lx.calc[,2]==i),1]	
	lx.temp <- lx.calc[(1+count):(6+count) ,1]
	lx.single[(i+2):(i+5)] <- coeffs.ages10_74 %*% lx.temp; ## matrix multiplication;	
	count <- count + 1
}	


#*** Gompertz parameters for ages 75+ ***

param.b <- ((log(lx.calc[18,1]) - log(lx.calc[17,1])) / (log(lx.calc[17,1]) - log(lx.calc[16,1])))^(1/5)
param.a <- exp((log(lx.calc[17,1]) - log(lx.calc[16,1])) / (param.b^75 * (param.b^5 - 1)))
param.c <- lx.calc[16,1] * exp((-param.b^75) * log(param.a))

remaining.ages <- 75:96
lx.single[76:length(lx.single)] <- param.c * (param.a^(param.b^remaining.ages))


#*** convert lx to Mx ***

dx <- lx.single[1:(length(lx.single)-1)] - lx.single[2:length(lx.single)]

dx <- c(dx,(lx.single[length(lx.single)]))

qx <- dx / lx.single

ax <- rep(0.5,97)
ax[1] <- 0.1

mx <- qx / (1-qx*(1-ax))


zzz <- data.frame(mx)
write.csv(zzz, file="d:/zzz.csv")
