# Impact of hormone therapy (HT) in the urinary microbiome of peri- or menopausal women
Final degree project for the Master Programme in Bioinformatics at Uppsala University. This project was carried out with the Research Group Johan Ärnlöv in the Department of Neurobiology, Care Sciences and Society of Karolinska Institute, under supervision of Tiscar Graells.
The data used for this analysis was obtained from the SCAPIS database, but following the GDPR rules, it will not be shared in this repository.

## Description and background of the project

Until recently, the urinary tract was considered sterile, and the presence of bacteria was interpreted as an indication of infection. However, advances in microbiological and sequencing techniques have demonstrated that the urinary tract harbors a diverse microbial community, whose composition may be associated with both urinary health and disease.

Urinary tract infections (UTIs) are the most common infections in women, and their risk increases significantly during and after menopause, often becoming recurrent (rUTIs). Women affected by these conditions are increasingly treated with oestrogen/progesterone hormone therapy, which appears to reduce the incidence of such infections.Nevertheless, the underlying biological mechanisms, and the ways in which hormone therapy may shape the urinary microbiome, remain largely unclear. It has been hypothesized that, similarly to the vaginal microbiome, hormone therapy may promote an increase in *Lactobacillus* species and a reduction in overall microbial diversity, resulting in a microbial profile similar to that of a healthy premenopausal woman.

For this project, bacterial abundance profiles from urine samples of peri- and postmenopausal women receiving hormone therapy were analyzed alongside matched control subjects. The aim was to assess and characterize the impact of hormone therapy (HT) on the composition and diversity of the urinary microbiome, using a case–control study design and statistical modeling approaches to identify potential treatment-associated differences.

## Workflow

The bacterial abundance of samples were given in relative abundance and centered-log ratio (CLR) transformed. In order to enable a more detailed assessment of the relationship between HT and the urinary microbiome,  the data were analyzed at three taxonomic levels: species, genus, and family.

As a preliminary step, diversity metrics were evaluated, including both alpha and beta diversity analyses. Following this, two linear regression models were fitted: a minimally adjusted model accounting only for batch effects, and a fully adjusted model including additional covariates to control for potential confounding.
The effect of HT duration was then examined by comparing continuous use (>4 months) versus recent use (<4 months), applying the same modeling framework.

All analyses were conducted in R v4.5.2 (R Foundation for Statistical Computing, Vienna, Austria) within the RStudio environment. 
