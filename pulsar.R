# load the data from file
pulsar <- read.table("Coding/R/pulsar-timing.txt", header=TRUE); #A table from J H Taylor and J M Weisberg, 1982, 
#describing pulsar timings. The paper has more details.
date <- pulsar$date; #The date each observation occurred.
timingVariation <- pulsar$dt; #How different this timing is from what would be expected if the orbital period was constant
error <- pulsar$error; #The standard deviation on the timing measurement.

#Define other variables
firstPulse = date[1]; #date of the first pulse
timeSinceFirstPulse <- (date - firstPulse); #time between each pulse and the first


#Define the first model. This model assumes there is no variation in timing (time between subsequent pulses is constant)
model.constant <- function(constants, timeSinceFirstPulse)
{
  timingVariation.mod <- constants[1]
  return(timingVariation.mod)
}

# Define the linear model. In this model, actual timing between pulses is constant, but the observed period
# is slightly different from the expected period, so timingDifference increases linearly with time
model.linear <- function(constants, timeSinceFirstPulse) 
{
  timingVariation.mod <- constants[1] + constants[2]*timeSinceFirstPulse
  return(timingVariation.mod)
}

# Define the quadratic model. This assume that the timing between pulses decreases over time due to gravitational interactions.
model.quadratic <- function(constants, timeSinceFirstPulse) 
{
  timingVariation.mod <- constants[1] + constants[2]*timeSinceFirstPulse + constants[3]*(timeSinceFirstPulse^2)
  return(timingVariation.mod)
}

#Define the likelihood model. This calculates the probability that the data would be observed if the model is true.
LogLikelihood <- function(constants, timeSinceFirstPulse, timingVariation, error, model) 
{
  timingVariation.mod <- model(constants, timeSinceFirstPulse)
  l <- dnorm(timingVariation, mean=timingVariation.mod, sd=error, log=TRUE) #This creates a probability density function with a normal distribution, one for each parameter. log=true means it outputs log(p)
  mlogl <- -sum(l) #log likelihood function is defined as the sum of the logs of the probability densities. We use -sum(l) so we can minimise rather than maximise, it's easier to calculate
  return(mlogl)
}

#Fitting the models
#Define the starting parameters. It doesn't particularly matter what these are, they're guesses based on the observed data and will be adjusted by R. Tested by altering slightly to make sure different values don't produce different results
constants.0 <- c(-1, -1, -1);
#Fit the constant model
result.constant <- optim(fn=LogLikelihood, method='Brent', par=constants.0[1],
                      hessian=TRUE, timeSinceFirstPulse=timeSinceFirstPulse, 
                      timingVariation=timingVariation, error=error, lower=constants.0[1]-1, upper=constants.0[1] + 1,
                      model=model.constant);
#Fit the linear model
result.linear <- optim(fn=LogLikelihood, par=c(constants.0[1], constants.0[2]),
                    hessian=TRUE, timeSinceFirstPulse=timeSinceFirstPulse, 
                    timingVariation=timingVariation, error=error,
                    model=model.linear);
#Fit the quadratic model
result.quadratic <- optim(fn=LogLikelihood, par=constants.0,
                     hessian=TRUE, timeSinceFirstPulse=timeSinceFirstPulse, 
                     timingVariation=timingVariation, error=error,
                     model=model.quadratic);
timingVariation.mod <- model.constant(result.constant$par, timeSinceFirstPulse);
timingVariation.mod2 <- model.linear(result.linear$par, timeSinceFirstPulse);
timingVariation.mod3 <- model.quadratic(result.quadratic$par, timeSinceFirstPulse);

#Plot everything
#Pick the colours for our 5 lines
colourData ="red";
colourNoVariation = "black";
colourConstant = "blue";
colourLinear = "darkgreen";
colourQuadratic = "blueviolet";
plot(date, timingVariation, col=colourData, type='o', xlab= "Date", 
     ylab = "Pulsar Timing Variation From Expected", main="Observed Pulsar Timings Plotted Against 3 Models");
segments(date, timingVariation-error, x1=date, y1=timingVariation+error, col=colourData);
abline(a=0, b=0, col= colourNoVariation); #This is what we would see if the pulsar had a constant period
abline(a=timingVariation.mod, b=0, col=colourConstant);
lines(date, timingVariation.mod2, col= colourLinear);
lines(date, timingVariation.mod3, col= colourQuadratic);
#Add a legend
legend('topright', legend=c("Data", "No Timing Variation", "Constant Model", "Linear Model",
                            "Quadratic Model"), lty=c(1,1),
        col=c(colourData, colourNoVariation, colourConstant, colourLinear, colourQuadratic));

#Next we find the p value
#Define Chi squared statistic
ChiSq <- function(constants, x, timingVariation, error, model) 
{
  timingVariation.mod <- model(constants, x) #Compute the y value for the given model
  x <- sum( ((timingVariation - timingVariation.mod)^2 )/ error^2) #For each point on a model we subtract it from an observed data point, square it, and divided by the standard deviation
  return(x)
}

#find chi squared for each model, using the parameters found earlier
x.mod1 <- ChiSq(result.constant$par, x, timingVariation, error,
                model=model.constant);
x.mod2 <- ChiSq(result.linear$par, x, timingVariation, error,
                model=model.linear);
x.mod3 <- ChiSq(result.quadratic$par, x, timingVariation, error,
                model=model.quadratic);

#N is simply the number of independent data points
N <- length(timingVariation);
degreesFreedom1 <- N - length(result.constant$par); #degrees of freedom for each model is found by subtracting the number of adjustable parameters from N
degreesFreedom2 <- N - length(result.linear$par);
degreesFreedom3 <- N - length(result.quadratic$par); #these are just the numbers 1, 2, 3

#pchisq is another built-in function from R, which computes a p value from the Chi squared statistic
p1 <- pchisq(x.mod1, df=degreesFreedom1, lower.tail=FALSE);
p2 <- pchisq(x.mod2, df=degreesFreedom2, lower.tail=FALSE);
p3 <- pchisq(x.mod3, df=degreesFreedom3, lower.tail=FALSE);

#Make a second estimate for the p value using simulations
xSequence <- seq(0, 3*degreesFreedom3, by=0.1); #creates a sequence of x values--these are the possible
timingVariation.pdf <- dchisq(xSequence, df=degreesFreedom3); #finds chi squared for each of these points
timingVariation.mod <- model.quadratic(result.quadratic$par, x);
N.sim <- 1000; #Number of simulations to run. If you use a larger number you get a more consistent result but it takes longer to calculate 
constants.sim <- array(0, dim=c(N.sim, length(result.quadratic$par))); #create an array with N rows and space for 3 parameters for each simulation
chisq.sim <- array(0, dim=N.sim); #create an array for all the chi squared values
for (i in 1:N.sim) 
{
  timingVariation.sim <- rnorm(length(timingVariation.mod), mean=timingVariation.mod, sd=error) #generate random data with a normal distribution with the same error as our data
  result.sim <- optim(fn=ChiSq, par=constants.0,
                      x=x, timingVariation=timingVariation.sim, error=error,
                      model=model.quadratic, control=list(parscale=constants.0)) #find the chi square value for each of our data sets
  constants.sim[i,] <- result.sim$par #add the calculated constants for each simulation to an array
  chisq.sim[i] <- result.sim$value #add the chi squared values to an array
}

print(paste("The p value for the quadratic model is " , mean(chisq.sim > x.min3), " using the simulation method"));
#what fraction of simulated chi squared values are greater than our actual value? This is the p value
print(paste("The p value for the quadratic model is " , p3, " using the pchisq function"))
#"paste" just concatenates the strings
