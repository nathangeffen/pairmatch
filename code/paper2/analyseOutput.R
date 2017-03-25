inp = read.csv("output_20k_3yr.csv", FALSE)
library(xtable)

confinterval <- function(vals) {
  tryCatch(
    {
      x = t.test(vals)
      x$conf.int
    },
    error = function(e) {
      NA
    }
  )
}


confinterval_02_5 <- function(vals) {
  tryCatch(
    {
      confinterval(vals)[1]
    },
      error = function(e) {
        NA
      }
    )
}

confinterval_97_5 <- function(vals) {
  tryCatch(
    {
      confinterval(vals)[2]
    },
    error = function(e) {
      NA
    }
  )
}

timings = inp[grepl("TIMING",inp$V2),]
timings = timings[grepl("AFTER",timings$V3),]
mean_timings = aggregate(timings$V5, by=list(timings$V1),FUN=mean)
print("TIMING")
print(mean_timings)

calcPrevalenceForDate <- function(date_of_analysis) {
  prevalences = inp[inp$V2 == "ANALYSIS",]
  prevalences = prevalences[grep("PREVALENCE|MALE",prevalences$V3),]
  prevalences = prevalences[prevalences$V5==date_of_analysis,]
  mean_prevalences =  aggregate(prevalences$V6, by=list(prevalences$V1,prevalences$V3),
                              FUN=mean)
  sd_prevalences = aggregate(prevalences$V6, by=list(prevalences$V1,prevalences$V3),
                           FUN=sd)
  ci_prevalences_02_5 = aggregate(prevalences$V6, by=list(prevalences$V1,prevalences$V3),
                           FUN=confinterval_02_5)
  ci_prevalences_97_5 = aggregate(prevalences$V6, by=list(prevalences$V1,prevalences$V3),
                               FUN=confinterval_97_5)
  output = data.frame(simulation=mean_prevalences$Group.1,
                    measure=mean_prevalences$Group.2, 
                    mean=mean_prevalences$x,
                    sd=sd_prevalences$x,
                    ci_02_5=ci_prevalences_02_5$x,
                    ci_97_5=ci_prevalences_97_5$x)
  output[is.na(output$mean)==FALSE,]

  #output[is.nan(output$mean)==FALSE,]
}

initialprevs = calcPrevalenceForDate(2017)
finalprevs = calcPrevalenceForDate(2020)

categories = c("ALL", "Males", "Females", "MSM", "WSW",
               "Male 15-19", "Female 15-19", "Male 20-24", "Female 20-24", 
               "Male 25-29", "Female 25-29", "Male 30-34", "Female 30-34", 
               "Male 35-39", "Female 35-39", "Male 40-44", "Female 40-44",
               "Male 45-49", "Female 45-49")
vals = rep(-1, length(categories))
tex_col_names = c("RPM 0.001", "RKPM 0.001", "CSPM 0.001",
                  "RPM 0.01", "RKPM 0.01", "CSPM 0.01",
                  "RPM 0.1", "RKPM 0.1", "CSPM 0.1")
r_col_names = gsub(" ","_", tex_col_names)
r_col_names = gsub("[.]", "_", r_col_names)
output_table = data.frame(row.names = categories)

tex_row_names = row.names(output_table)
s = toupper(tex_row_names)
s = gsub(" ","_0", s)
s = gsub("-","-0", s)
s = gsub("MALE", "MALE_AGE", s)
s = gsub("WSW", "WSWPREVALENCE", s)
s = gsub("MSM", "MSMPREVALENCE", s)
s = gsub("ALL", "PREVALENCE", s)
r_row_names = gsub("MALE_AGES", "MALEPREVALENCE", s)
for (i in c(0:length(r_row_names))) {
  for (j in c(0:length(r_col_names))) {
    expected = format(round(subset(finalprevs,measure==r_row_names[i] & 
                            simulation==r_col_names[j])$mean[1], 3))
    stddev = format(round(subset(finalprevs,measure==r_row_names[i] & 
                              simulation==r_col_names[j])$sd[1], 3))
    c_025 = format(round(subset(finalprevs,measure==r_row_names[i] & 
                                   simulation==r_col_names[j])$ci_02_5[1], 3)) 
    c_975 = format(round(subset(finalprevs,measure==r_row_names[i] & 
                                  simulation==r_col_names[j])$ci_97_5[1], 3))      
    output_table[tex_row_names[i],  tex_col_names[j] ] = paste(expected, " ", 
                                                               "[", c_025, ";", c_975,"]",
                                                               sep="") 
  }
}

xtable(output_table)
# cspm_0_01 = prevalences[prevalences$V1=="CSPM_0_01",]
# cspm_0_01 = cspm_0_01[cspm_0_01$V3=="PREVALENCE",]
# hist(cspm_0_01$V6)
# 
# rpm_0_01 = prevalences[prevalences$V1=="RPM_0_01",]
# rpm_0_01 = rpm_0_01[rpm_0_01$V3=="PREVALENCE",]
# hist(rpm_0_01$V6)
# 
# rpm_0_1 = prevalences[prevalences$V1=="RPM_0_1",]
# rpm_0_1 = rpm_0_1[rpm_0_1$V3=="PREVALENCE",]
# hist(rpm_0_1$V6)