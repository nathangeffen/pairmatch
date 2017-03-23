inp = read.csv("output2.csv", FALSE)
timings = inp[grepl("TIMING",inp$V2),]
timings = timings[grepl("AFTER",timings$V3),]
mean_timings = aggregate(timings$V5, by=list(timings$V1),FUN=mean)
print("TIMING")
print(mean_timings)

prevalences = inp[inp$V3 == "PREVALENCE", ]
prevalences = prevalences[prevalences$V5==2020.000,]
mean_prevalences = aggregate(prevalences$V6, by=list(prevalences$V1),
                             FUN=mean)
print("OVERALL PREVALENCE")
print(mean_prevalences)

prevalences = inp[inp$V3 == "MALEPREVALENCE", ]
prevalences = prevalences[prevalences$V5==2020.000,]
mean_prevalences = aggregate(prevalences$V6, by=list(prevalences$V1),
                             FUN=mean)
print("MALE PREVALENCE")
print(mean_prevalences)

prevalences = inp[inp$V3 == "FEMALEPREVALENCE", ]
prevalences = prevalences[prevalences$V5==2020.000,]
mean_prevalences = aggregate(prevalences$V6, by=list(prevalences$V1),
                             FUN=mean)
print("FEMALE PREVALENCE")
print(mean_prevalences)

prevalences = inp[inp$V3 == "MSMPREVALENCE", ]
prevalences = prevalences[prevalences$V5==2020.000,]
mean_prevalences = aggregate(prevalences$V6, by=list(prevalences$V1),
                             FUN=mean)
print("MSM PREVALENCE")
print(mean_prevalences)

prevalences = inp[inp$V3 == "WSWPREVALENCE", ]
prevalences = prevalences[prevalences$V5==2020.000,]
mean_prevalences = aggregate(prevalences$V6, by=list(prevalences$V1),
                             FUN=mean)
print("WSW PREVALENCE")
print(mean_prevalences)
