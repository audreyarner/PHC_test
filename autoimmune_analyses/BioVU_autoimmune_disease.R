setwd("~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/RotationProject/")
Codes = read.delim("auto_immune_codes_for_people_with_labor_recs.txt")

#remove individuals with juvenile diseases
Adult = Codes[!grepl("juvenile|Juvenile", Codes$CONCEPT_NAME), ]

Preterm = read.delim("preterm_births_for_all_people.txt")
Term = read.delim("term_births_for_all_people.txt")
Abortion = read.delim("abortions_or_miscarriages_for_all_people.txt")

#remove duplicated individuals
UniquePreterm = as.data.frame(Preterm[!duplicated(Preterm$PERSON_SOURCE_VALUE),])
UniqueTerm = as.data.frame(Term[!duplicated(Term$PERSON_SOURCE_VALUE),])
UniqueAbortion = as.data.frame(Abortion[!duplicated(Abortion$PERSON_SOURCE_VALUE),])

#merge data
Together = merge(UniqueTerm, UniquePreterm, by=c("PERSON_SOURCE_VALUE", "PERSON_ID", "FORM_NAME"))
Together = merge(Together, Abortion, by=c("PERSON_SOURCE_VALUE", "PERSON_ID", "FORM_NAME"))

#in the following lines, make new rows with number of times each individual had each type of pregnancy
Together$Term = 0
Together$Preterm = 0
Together$Abortion = 0
Together$Autoimmune = 0
Together$KidAutoimmune = 0

for(i in 1:nrow(Together)){
  if(Together$FIELD_VALUE.x[i] %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
    Together$Term[i] = Together$FIELD_VALUE.x[i]
  }
  if(Together$FIELD_VALUE.y[i] %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")){
    Together$Preterm[i] = Together$FIELD_VALUE.y[i]
  }
  if(Together$FIELD_VALUE[i] %in% c("0","1", "2", "3", "4", "5", "6", "7", "8", "9")){
    Together$Abortion[i] = Together$FIELD_VALUE[i]
  }
  if(Together$PERSON_SOURCE_VALUE[i] %in% Adult$PERSON_SOURCE_VALUE){
    Together$Autoimmune[i] = 1
  }
  if(Together$PERSON_SOURCE_VALUE[i] %in% Codes$PERSON_SOURCE_VALUE){
    Together$KidAutoimmune[i] = 1
  }
}


#make everything numeric
Together$Term = as.numeric(Together$Term)
Together$Preterm = as.numeric(Together$Preterm)
Together$Abortion = as.numeric(Together$Abortion)

#two ways to count pregnancies
Together$TotalPreg = Together$Term + Together$Preterm
Together$AllPreg = Together$Term + Together$Preterm + 0.5*Together$Abortion

#make column for binary of if had at least one pregnancy or not
for(i in 1:nrow(Together)){
  if(Together$TotalPreg[i]!= 0){
    Together$BinaryPreg[i] = 1
  } else(Together$BinaryPreg[i]=0)
}

#read in demographic information
Demographic = read.table("demographic_info_for_people_with_autoimmune_or_labor_recs.txt")
#find age, given the data was pulled in 2023
Demographic$age = 2023-as.numeric(gsub("\\-.*","",Demographic$PERSON_ID))
DemoAlso = merge(Together, Demographic, by.x = "PERSON_ID", by.y = "PERSON_SOURCE_VALUE")

#models to determine whether autoimmune disease is predicted by pregnancy
Model1 = glm(Autoimmune~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model1)$coefficients[2]
summary(Model1)$coefficients[2,4] 

#test whether juvenile autoimmune diseases impact results
Juvenile = glm(KidAutoimmune~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Juvenile)$coefficients[2]
summary(Juvenile)$coefficients[2,4] 

#only include pregnant individuals
OnlyPreg = subset(DemoAlso, TotalPreg !=0)
Model2 = glm(Autoimmune~age + ethnicity + EHR_reported_race, data = OnlyPreg, family = binomial)
summary(Model2)$coefficients[2]
summary(Model2)$coefficients[2,4] 

#model with binary of ever pregnant or not as the predictor
Binary = glm(TotalPreg~Autoimmune + age + ethnicity + EHR_reported_race, data = DemoAlso, family = gaussian)
summary(Binary)$coefficients[2]
summary(Binary)$coefficients[2,4] 



CodeTable = as.data.frame(table(Adult$CONCEPT_NAME))
write.table(CodeTable, "AutoimmuneICDFreq.txt", sep = "\t", quote = F, row.names = F)
length(unique(Adult$CONCEPT_NAME))
length(unique(Adult$PERSON_SOURCE_VALUE))

lupus = Adult[grep("lupus", Adult$CONCEPT_NAME), ]
lupusSmall = as.data.frame(unique(lupus[c("CONCEPT_NAME", "CONCEPT_CODE")]))
lupusSmall$NumberIndivids = length(unique(lupus$PERSON_SOURCE_VALUE))
lupusSmall$UmbrellaTerm = "SystemicLupusErythematosus"


heumatoid = Adult[grep("heumatoid", Adult$CONCEPT_NAME), ]
heumatoidSmall = as.data.frame(unique(heumatoid[c("CONCEPT_NAME", "CONCEPT_CODE")]))
heumatoidSmall$NumberIndivids = length(unique(heumatoid$PERSON_SOURCE_VALUE))
heumatoidSmall$UmbrellaTerm = "RheumatoidArthritis"

Goiter = Adult[grep("goiter", Adult$CONCEPT_NAME), ]
GoiterSmall = as.data.frame(unique(Goiter[c("CONCEPT_NAME", "CONCEPT_CODE")]))
GoiterSmall$NumberIndivids = length(unique(Goiter$PERSON_SOURCE_VALUE))
GoiterSmall$UmbrellaTerm = "GravesDisease"

colitis = Adult[grep("colitis|proctitus|Ulcerative", Adult$CONCEPT_NAME), ]
colitisSmall = as.data.frame(unique(colitis[c("CONCEPT_NAME", "CONCEPT_CODE")]))
colitisSmall$NumberIndivids = length(unique(colitis$PERSON_SOURCE_VALUE))
colitisSmall$UmbrellaTerm = "UlcerativeColitis"

soriasis = Adult[grep("soriasis|psoriatic", Adult$CONCEPT_NAME), ]
soriasisSmall = as.data.frame(unique(soriasis[c("CONCEPT_NAME", "CONCEPT_CODE")]))
soriasisSmall$NumberIndivids = length(unique(soriasis$PERSON_SOURCE_VALUE))
soriasisSmall$UmbrellaTerm = "Psoriasis"

cirrhosis = Adult[grep("cirrhosis", Adult$CONCEPT_NAME), ]
cirrhosisSmall = as.data.frame(unique(cirrhosis[c("CONCEPT_NAME", "CONCEPT_CODE")]))
cirrhosisSmall$NumberIndivids = length(unique(cirrhosis$PERSON_SOURCE_VALUE))
cirrhosisSmall$UmbrellaTerm = "BiliaryCirrhosis"

thyroiditis = Adult[grep("thyroiditis", Adult$CONCEPT_NAME), ]
thyroiditisSmall = as.data.frame(unique(thyroiditis[c("CONCEPT_NAME", "CONCEPT_CODE")]))
thyroiditisSmall$NumberIndivids = length(unique(thyroiditis$PERSON_SOURCE_VALUE))
thyroiditisSmall$UmbrellaTerm = "HashimotoThyroiditis"

Celiac = Adult[grep("Celiac", Adult$CONCEPT_NAME), ]
CeliacSmall = as.data.frame(unique(Celiac[c("CONCEPT_NAME", "CONCEPT_CODE")]))
CeliacSmall$NumberIndivids = length(unique(Celiac$PERSON_SOURCE_VALUE))
CeliacSmall$UmbrellaTerm = "CeliacDisease"

sclerosis = Adult[grep("sclerosis", Adult$CONCEPT_NAME), ]
sclerosisSmall = as.data.frame(unique(sclerosis[c("CONCEPT_NAME", "CONCEPT_CODE")]))
sclerosisSmall$NumberIndivids = length(unique(sclerosis$PERSON_SOURCE_VALUE))
sclerosisSmall$UmbrellaTerm = "MultipleSclerosis"

Crohns = Adult[grep("Crohn's|enteritis", Adult$CONCEPT_NAME), ]
CrohnsSmall = as.data.frame(unique(Crohns[c("CONCEPT_NAME", "CONCEPT_CODE")]))
CrohnsSmall$NumberIndivids = length(unique(Crohns$PERSON_SOURCE_VALUE))
CrohnsSmall$UmbrellaTerm = "CrohnsDisease"

spondylitis = Adult[grep("spondylitis", Adult$CONCEPT_NAME), ]
spondylitisSmall = as.data.frame(unique(spondylitis[c("CONCEPT_NAME", "CONCEPT_CODE")]))
spondylitisSmall$NumberIndivids = length(unique(spondylitis$PERSON_SOURCE_VALUE))
spondylitisSmall$UmbrellaTerm = "AnkylosingSpondylitis"

gravis = Adult[grep("gravis", Adult$CONCEPT_NAME), ]
gravisSmall = as.data.frame(unique(gravis[c("CONCEPT_NAME", "CONCEPT_CODE")]))
gravisSmall$NumberIndivids = length(unique(gravis$PERSON_SOURCE_VALUE))
gravisSmall$UmbrellaTerm = "MyastheniaGravis"

diabetes = Adult[grep("diabetes", Adult$CONCEPT_NAME), ]
diabetesSmall = as.data.frame(unique(diabetes[c("CONCEPT_NAME", "CONCEPT_CODE")]))
diabetesSmall$NumberIndivids = length(unique(diabetes$PERSON_SOURCE_VALUE))
diabetesSmall$UmbrellaTerm = "Type1Diabetes"

AllAutoimmune = rbind(diabetesSmall, gravisSmall, spondylitisSmall, CrohnsSmall,
                      sclerosisSmall, CeliacSmall, thyroiditisSmall, cirrhosisSmall,
                      soriasisSmall, colitisSmall, GoiterSmall, lupusSmall, heumatoidSmall)

AutoimmuneLarge = as.data.frame(unique(AllAutoimmune[c("UmbrellaTerm", "NumberIndivids")]))

write.table(AllAutoimmune, "BioVUAutoimmuneUmbrella.txt", quote = F, row.names=F, sep = "\t")

#make new rows for whether each individual has had each autoimmune disease
DemoAlso$diabetes = 0
DemoAlso$gravis = 0
DemoAlso$spondylitis = 0
DemoAlso$crohns = 0
DemoAlso$sclerosis = 0
DemoAlso$Celiac = 0
DemoAlso$thyroiditis = 0
DemoAlso$cirrhosis = 0
DemoAlso$psoriasis = 0
DemoAlso$colitis = 0
DemoAlso$Goiter = 0
DemoAlso$Lupus = 0
DemoAlso$heumatoid = 0

for(i in 1:nrow(DemoAlso)){
  if(DemoAlso$PERSON_SOURCE_VALUE[i] %in% colitis$PERSON_SOURCE_VALUE){
    DemoAlso$colitis[i] = 1
  }
  if(DemoAlso$PERSON_SOURCE_VALUE[i] %in% soriasis$PERSON_SOURCE_VALUE){
    DemoAlso$psoriasis[i] = 1
  }
  if(DemoAlso$PERSON_SOURCE_VALUE[i] %in% Crohns$PERSON_SOURCE_VALUE){
    DemoAlso$crohns[i] = 1
  }
  if(DemoAlso$PERSON_SOURCE_VALUE[i] %in% gravis$PERSON_SOURCE_VALUE){
    DemoAlso$gravis[i] = 1
  }
  if(DemoAlso$PERSON_SOURCE_VALUE[i] %in% diabetes$PERSON_SOURCE_VALUE){
    DemoAlso$diabetes[i] = 1
  }
  if(DemoAlso$PERSON_SOURCE_VALUE[i] %in% lupus$PERSON_SOURCE_VALUE){
    DemoAlso$Lupus[i] = 1
  }
  if(DemoAlso$PERSON_SOURCE_VALUE[i] %in% thyroiditis$PERSON_SOURCE_VALUE){
    DemoAlso$thyroiditis[i] = 1
  }
  if(DemoAlso$PERSON_SOURCE_VALUE[i] %in% Goiter$PERSON_SOURCE_VALUE){
    DemoAlso$Goiter[i] = 1
  }
  if(DemoAlso$PERSON_SOURCE_VALUE[i] %in% sclerosis$PERSON_SOURCE_VALUE){
    DemoAlso$sclerosis[i] = 1
  }
  if(DemoAlso$PERSON_SOURCE_VALUE[i] %in% heumatoid$PERSON_SOURCE_VALUE){
    DemoAlso$heumatoid[i] = 1
  }
  if(DemoAlso$PERSON_SOURCE_VALUE[i] %in% cirrhosis$PERSON_SOURCE_VALUE){
    DemoAlso$cirrhosis[i] = 1
  }
  if(DemoAlso$PERSON_SOURCE_VALUE[i] %in% Celiac$PERSON_SOURCE_VALUE){
    DemoAlso$Celiac[i] = 1
  }
  if(DemoAlso$PERSON_SOURCE_VALUE[i] %in% spondylitis$PERSON_SOURCE_VALUE){
    DemoAlso$spondylitis[i] = 1
  }
}

#break models up by specific disease
Model3 = glm(diabetes~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model3)$coefficients[2,4]

Model3 = glm(gravis~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model3)$coefficients[2,4]

Model3 = glm(crohns~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model3)$coefficients[2,4]

Model3 = glm(psoriasis~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model3)$coefficients[2,4]

Model3 = glm(colitis~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model3)$coefficients[2,4]

Model3 = glm(Lupus~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model3)$coefficients[2,4]

Model4 = glm(thyroiditis~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model4)$coefficients[2,4]

Model5 = glm(Goiter~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model5)$coefficients[2,4]

Model6 = glm(sclerosis~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model6)$coefficients[2,4]

Model7 = glm(heumatoid~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model7)$coefficients[2,4]

Model8 = glm(cirrhosis~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model8)$coefficients[2,4]

Model9 = glm(Celiac~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model9)$coefficients[2,4]

Model9 = glm(spondylitis~TotalPreg + age + ethnicity + EHR_reported_race, data = DemoAlso, family = binomial)
summary(Model9)$coefficients[2,4]

