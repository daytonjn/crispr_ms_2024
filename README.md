# CRISPR manuscript

**Description:** Scripts to reproduce analyses reported in Dayton et al. ([2024](https://doi.org/10.1111/imb.12959)): Efficient CRISPR/Cas9-mediated genome editing in the European corn borer, *Ostrinia nubilalis*.

**Purpose:** A CRISPR/Cas9 gene-editing & microinjection protocol was developed for the European corn borer, a major agricultural pest in the USA. The contribution of *period*, *pigment-dispersing factor receptor*, and genetic background on daily behaviors (eclosion timing) & development (diapause propensity, mass) were quantified. Analyses were all run in R.  

**Analysis Workflow:**
 1) Quantified editing efficiencies & outcomes of different CRISPR/Cas9 gene-editing protocols (low vs. high Cas9/sgRNA concentration)
    - *editing_statistics.R*
 2) Evaluated effect of *per*/*pdfr* mutagenesis in injected-F<sub>0</sub> vs. wild-type F<sub>0</sub> on rhythmic strength (Winfree's R) & distribution of daily behavior timing (eclosion)
    - *figure2_01.Rmd* & *figure2_02.R*
    - *figureS1_eclosion.Rmd* & *figureS2_eclosion.Rmd*
 3) Leveraged Quantitative Complementation Test (QCT) approach and mutagenesis to evaluate evidence for additive effects of *ptth* knockout, *per* genotype (uni- vs. bivoltine), and a **genetic** interaction (epistasis) between *ptth* AND *per* on diapause incidence & mass (see below). 
    - *figure4_QCT_model_schematic.R*
    - *figure5_QCT_experiment.Rmd*

**Manuscript:** *dayton_2024_ms.pdf*
