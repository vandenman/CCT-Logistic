### General

- [ ] Probably rename machine_learning.R to something else.
- [ ] Move functions in machine_learning.R to separate file.
- [ ] Look into more elaborate performance metrics.
- [ ] Look into multiple training/ test sets?
- [ ] There are 20 patients lost due to missing values somewhere. This is a 20% loss, a bit much.
- [ ] Look into predictions based on partial information (e.g., only patient characteristics, no IFBE).


### Machine Learning

- [x] Impute missing values?
    - [x] used mice for imputations


### LTRM

- [x] Run normal LTRM on the items. Somehow check if results make sense.
- [ ] Add logistic regression to LTRM to do predictions.
- [ ] Find references on IFBE questionnaire.
- [x] try the QR trick for the sum to zero constraint rather than what I do now

- [ ] S3 method for converting to Stan.
- [ ] add implied thresholds instead of "free_thresholds".


### Discuss with EJ

- [ ] Approach where the 4 current methods use "stupidly" average the ratings
- [ ] Approaches for predictions based on partial information


### TODO

- [x] combine logistic with non-logistic model into one Stan file?
- [ ] two time points!
- [ ] additional patient parameter

- [ ] add bool for running logistic regression or not

