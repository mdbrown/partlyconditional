
load("/Volumes/zheng_y/PASS/Anna/PSA Kinetics/PASS clinical models/For Marshall_risk calculator/alldata_nomisspsa")

alldata_nomisspsa$inter<-ifelse(alldata_nomisspsa$fivealpha_ever==0,
                                0,alldata_nomisspsa$time_since_dx)

alldata_nomisspsa$inter<-ifelse(alldata_nomisspsa$fivealpha_ever==0,
                                0,alldata_nomisspsa$time_since_dx)

mylme <- lme(log_psa~1+time_since_dx+fivealpha_ever+inter,
             data = alldata_nomisspsa,
             random = ~1+time_since_dx|study_participant_id)

myblup <- BLUP(marker = "log_psa",
               measurement.time = "time_since_dx",
               fixed = c("time_since_dx", "fivealpha_ever", "inter"),
               random = c("time_since_dx"),
               id = "study_participant_id" ,
               data = alldata_nomisspsa )


newdata <- alldata_nomisspsa %>%
  filter(is.element(study_participant_id, c("8800906", "8800914")))


predict(myblup, newdata = newdata)



