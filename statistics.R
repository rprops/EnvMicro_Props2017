### Generalized mixed model to compare alpha-div and cell density for treatment effect
library("lme4")
library("lmtest")
test.data <- aggregate(D2~Sample/Treatment/Time, data=results, FUN = mean)

m2 <- lmer(D2 ~ Treatment + (Time|Sample), test.data)
plot(m2)

test2 <- lme(fixed=D2~Treatment , random=~1|Time/Sample, data=test.data, correlation=corAR1(form=~Time|Sample))
plot(test2)

# test <- lme(fixed=D1~Treatment, random=~1|bio_rep, data=results)
# anova(test)
# shapiro.test(residuals(test))
# bptest(residuals(test)~fitted(test))
# plot(residuals(test)~fitted(test))
# qqnorm(residuals(test))
# qqline(residuals(test))
# 
# test <- lme(fixed=D0~Treatment, random=~1|bio_rep, data=results)
# anova(test)
# shapiro.test(residuals(test))
# bptest(residuals(test)~fitted(test))
# plot(residuals(test)~fitted(test))
# qqnorm(residuals(test))
# qqline(residuals(test))

## Cell density

test <- lme(fixed=Total.cells~Treatment, random=~1|bio_rep, data=results)
anova(test)
shapiro.test(residuals(test))
bptest(residuals(test)~fitted(test))
plot(residuals(test)~fitted(test))
qqnorm(residuals(test))
qqline(residuals(test))

### Test for effect between two different treatments
