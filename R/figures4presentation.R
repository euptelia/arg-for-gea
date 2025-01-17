# Plots for presentation

# Recurrent environmental change/range change
# cycle = seq(100000)
cycle = seq(20000)
envPeriod = 10000
peakTime=0.2*envPeriod
expansionSpeed=0.9/peakTime
contractionSpeed=0.9/(envPeriod-peakTime)
currentPhase = cycle %% envPeriod

b = ifelse(currentPhase < peakTime,
       0.1+currentPhase*expansionSpeed,
       1-(currentPhase-peakTime)*contractionSpeed)
plot(cycle, b,
     lab = c(10, 8, 7),
     xlab = "Generations",
     ylab = "Upper boundary of x")