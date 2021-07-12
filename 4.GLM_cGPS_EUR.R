#######################################################################
##################### 4. GLM_cognitive gps #####################
#######################################################################

##################### visualization default setting ##################### 
m.col = "#5B84B1FF"
f.col = "#FC766AFF"
sex.col = scale_fill_discrete(type = c(m.col, f.col))
sex.coll = scale_color_discrete(type = c(m.col, f.col))
sex.scale = scale_x_discrete(labels = c("Male", "Female"))
my.theme <- theme_bw() + theme(legend.position = "none", axis.title = element_text(size = 20), 
                               axis.text = element_text(size = 15, color = "black"))
gp <- ggplot(data = df) + sex.col

##################### load data #####################
########### total ########### 
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/connectome21/2021-1/project/Sex/analysis data")
gps.total <- read.csv("demo.gps.r.total.csv", na = "")
gps.total$sex <- as.factor(gps.total$sex)

gps.total$sex <- as.factor(gps.total$sex)
gps.total$race.ethnicity <- as.factor(gps.total$race.ethnicity)
gps.total$high.educ <- as.integer(gps.total$high.educ)
#gps.total$income <- as.factor(gps.total$income)
gps.total$married <- as.factor(gps.total$married)
gps.total$abcd_site <- as.factor(gps.total$abcd_site)

library(dplyr)
gps.total.1 <- gps.total %>% filter(sex ==1)
gps.total.2 <- gps.total %>% filter(sex == 2)

########### train ########### 
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/connectome21/2021-1/project/Sex/analysis data")
gps.train <- read.csv("demo.gps.r.train.csv", na = "")
gps.train$sex <- as.factor(gps.train$sex)

#gps.train <- mutate(gps.train, morctsex.1z = (morctsex.1 - mean(morctsex.1)) / sd(morctsex.1))
#gps.train <- mutate(gps.train, morctsex.1log = cuberoot(morctsex.1))

gps.train$sex <- as.factor(gps.train$sex)
gps.train$race.ethnicity <- as.factor(gps.train$race.ethnicity)
gps.train$high.educ <- as.integer(gps.train$high.educ)
#gps.train$income <- as.factor(gps.train$income)
gps.train$married <- as.factor(gps.train$married)
gps.train$abcd_site <- as.factor(gps.train$abcd_site)

library(dplyr)
gps.train.1 <- gps.train %>% filter(sex ==1)
gps.train.2 <- gps.train %>% filter(sex == 2)

##################### load family ID #####################
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/connectome21/2021-1/project/Sex/analysis data")
familyID <- readxl::read_xlsx("familyID.xlsx", na = "") #11878
familyID.uniq <- familyID[!duplicated(familyID[,"rel_family_id"]), ] #9856

##################### merge GPS and family ID #####################
gps.train.fm <- merge(familyID, gps.train, by = 'subjectkey')
gps.train.fm.1 <- merge(familyID, gps.train.1, by = 'subjectkey')
gps.train.fm.2 <- merge(familyID, gps.train.2, by = 'subjectkey')

#################################### 
############ EUR ############ 
#################################### 
###EUR only
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/abcd final/Gene/PCA")
EUR <- read.csv("ABCD_QCed_2021_PCair_8620samples_EUR.sample", sep = "\t")
EUR$FID <- gsub("AB000", "", EUR$FID)
EUR$FID <- gsub("[0-9][0-9][0-9][0-9]_", "", EUR$FID)
names(EUR)[1] <- 'subjectkey'

###df.eur
df.eur <- merge(EUR[1], gps.train.fm, by = 'subjectkey')
df.eur.1 <- merge(EUR[1], gps.train.fm.1, by = 'subjectkey')
df.eur.2 <- merge(EUR[1], gps.train.fm.2, by = 'subjectkey')

##################### GLM default setting#####################
library(ggExtra)
library(lmerTest)
library(MASS)

###### GLM of cGPS In Both Sex Group ######
###EA
df.eur.1$rel_family_id <- as.factor(df.eur.1$rel_family_id)
gps.train.EA <- lm(data = df.eur, 
                    formula = EA ~ morctsex.1 + sex + age + height + weight + BMI + high.educ + income + 
                      married + abcd_site)
summary(gps.train.EA)

EA.eta <- effectsize::eta_squared(aov(gps.train.EA, data = df.eur), partial = T, ci = 0.95)[1,2] #family=binomial
EA.eta

###CP
gps.train.CP <- lm(data = df.eur, 
                     formula = CP ~ morctsex.1 + sex + age + height + weight + BMI + high.educ + income + 
                       married + abcd_site + (1 | rel_family_id))
summary(gps.train.CP)
CP.eta <- effectsize::eta_squared(aov(gps.train.CP, data = df.eur), partial = T, ci = 0.95)[1,2] #family=binomial
CP.eta

###IQ
gps.train.IQ <- lm(data = df.eur, 
                    formula = IQ ~ morctsex.1 +  sex + age + height + weight + BMI + high.educ + income + 
                      married + abcd_site)
summary(gps.train.IQ)

IQ.eta <- effectsize::eta_squared(aov(gps.train.IQ, data = df.eur), partial = T, ci = 0.95)[1,2] #family=binomial
IQ.eta


###### GLM of cGPS In Males ######
###EA
gps.train.EA.1 <- lm(data = df.eur.1, 
                     formula = EA ~ morctsex.1 + age + height + weight + BMI + high.educ + income + 
                       married  + abcd_site)
summary(gps.train.EA.1)

EA.eta.1 <- effectsize::eta_squared(aov(gps.train.EA.1, data = df.eur.1), partial = T, ci = 0.95)[1,2] #family=binomial
EA.eta.1

###CP
gps.train.CP.1 <- lm(data = df.eur.1, 
                     formula = CP ~ morctsex.1 + age + height + weight + BMI + high.educ + income + 
                       married + abcd_site)
summary(gps.train.CP.1)
CP.eta.1 <- effectsize::eta_squared(aov(gps.train.CP.1, data = df.eur.1), partial = T, ci = 0.95)[1,2] #family=binomial
CP.eta.1

###IQ
gps.train.IQ.1 <- lm(data = df.eur.1, 
                     formula = IQ ~ morctsex.1 + age + height + weight + BMI + high.educ + income + 
                       married + abcd_site)
summary(gps.train.IQ.1)

IQ.eta.1 <- effectsize::eta_squared(aov(gps.train.IQ.1, data = df.eur.1), partial = T, ci = 0.95)[1,2] #family=binomial
IQ.eta.1

###### GLM of cGPS In Females ######
###EA
gps.train.EA.2 <- lm(data = df.eur.2, 
                       formula = EA ~ morctsex.1 + age + height + weight + BMI + high.educ + income + 
                         married + abcd_site)
summary(gps.train.EA.2)

EA.eta.2 <- effectsize::eta_squared(aov(gps.train.EA.2, data = df.eur.2), partial = T, ci = 0.95)[1,2] #family=binomial
EA.eta.2

###CP
gps.train.CP.2 <- lm(data = df.eur.2, 
                       formula = CP ~ morctsex.1 + age + height + weight + BMI + high.educ + income + 
                         married + abcd_site)
summary(gps.train.CP.2)

CP.eta.2 <- effectsize::eta_squared(aov(gps.train.CP.2, data = df.eur.2), partial = T, ci = 0.95)[1,2] #family=binomial
CP.eta.2

###IQ
gps.train.IQ.2 <- lm(data = df.eur.2, 
                       formula = IQ ~ morctsex.1 + age + height + weight + BMI + high.educ + income + 
                         married + abcd_site)
summary(gps.train.IQ.2)

IQ.eta.2 <- effectsize::eta_squared(aov(gps.train.IQ.2, data = df.eur.2), partial = T, ci = 0.95)[1,2] #family=binomial
IQ.eta.2

