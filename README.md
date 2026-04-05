# Spatial-IPM-Telemetry
Jolly-Seber spatial capture-recapture model with integrated telemetry survival data

This model adds a telemetry observation process to the Jolly-Seber model here (JS-SCR-Dcov):

https://github.com/benaug/Jolly-Seber-N-Prior-DA

The model is set up with 2 SCR observation processes, Method 1, which is used to deploy telemetry collars, and Method 2, 
which is another arbitrary capture or detection process that also provides individual ID. Methods 1 and 2 can be deployed
in different years. The test script is set up with 7 years where collars are deployed in years 2 and 5 and the other 
observation process is used in years 1,4, and 7. There is no sampling in years 3 and 6, but you may have telemetry data for
previously collared individuals that are still alive. Data simulator allows you to control how long telemetry collars last.

Telemetry survival data is assumed to be uninformatively censored. Telemetry only provides more survival information than the SCR
data to the extent individuals are known to be alive longer via the collar than via SCR. 

To include telemetry survival information without introducing bias into the survival parameter estimates, you need to model the marking
process to explain why dead individuals do not receive collars. The marking process does not need to be spatial unless survival 
varies spatially, which you may want to consider. However, this model includes a spatial density covariate and you need to model the
marking process to prevent bias in the spatial density covariate parameter estimates, unless you collar individuals at random with
respect to space (no one does), or you don't allow the telemetry data to inform their activity centers. Therefore, this model
assumes there is a spatial marking process and that data was recorded.