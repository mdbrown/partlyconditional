
Looking at `newd` further, we see that subject with id 9 only has marker data up to 18 months, and subject 74 has information up to time 12. If we wish to make predictions at measurement time 24 for these participants, one option is to 'carry forward' the information from the previous visit to a new measurement time.   We can use the helper function `carry_forward` to accomplish this.

```{r}
newd.2 <- carry_forward(newd,
                        id = "sub.id",
                        measurement.time = "meas.time",
                        new.times = c(24))
newd.2
```

Now we see that new rows are added for participant 9 and 74 with meas.time = 24. All other information (including predictors, status, and survival time) are carried forward from the last measurement time observed for each subject. *Unfortunately, this also includes the previous marker BLUP estimate that was calculated using the previous measurement time.* To update the blup estimates for each marker, we call `predict` on the blup model object (below it is called `myblup.marker1`). This step would need to be repeated for each blup marker model if more than one is used.


```{r}
#update blups for new measurement time by calling predict on the blup model
newd.2$marker_1_blup <- predict(myblup.marker1, newd.2)$fitted.blup
```

Note: The above step is only necessary if pc models employ BLUP-smoothed markers as predictors.
