GSL-rl : Gulf of St. Lawrence reference library
================

Reference sequences for metabarcording species assignments in the Gulf
of St. Lawrence

**Main maintainer:** Audrey Bourret  
**Affiliation:** Fisheries and Oceans Canada (DFO)  
**Group:** Laboratory of genomics  
**Location:** Maurice Lamontagne Institute  
**Affiliated publication:** Maximizing the reliability and the number of
species assignments in metabarcoding studies (in prepatation). Bourret,
A., Nozères, C., Parent, É., Parent, G.J.  
**Contact:** <audrey.bourret@dfo-mpo.gc.ca>

-   [Description of GSL-rl](#description-of-gsl-rl)
-   [Status](#status)
-   [Contents](#contents)
    -   [Subsections within contents](#subsections-within-contents)
-   [Requirements](#requirements)
-   [Caveats](#caveats)
-   [Uncertainty](#uncertainty)
-   [Acknowledgements](#acknowledgements)
-   [References](#references)

## Description of the GSL-rl

GSL-rl is a collection of curated and annotated sequences for performing
genetic assignments using COI barcodes for marine fauna in the Gulf of
St. Lawrence. It actually covered 651 targeted species, and sequences
are currently available for **439 species**.

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

## Status

“Ongoing-improvements”

## Contents

### Folder structure

    .
    ├── GSL-rl_COI     # Main folder containing reference sequences for COI 
    └── README.md

### GSL-rl COI

Two fasta files are available, one for the full fragment (\~650 pb;
Folmer) and another one cut to fit common smaller metabarcoding region
(\~313 pb; Leray). The smaller fragment is also composed of less
sequences because exact duplicated sequences were removed.

-   GSL-rl\_COI\_Folmer650pb\_450taxa\_1304seq.fasta  
-   GSL-rl\_COI\_Leray313pb\_450taxa\_1040seq.fasta

``` r
DNA.folmer <- Biostrings::readDNAStringSet("GSL-rl_COI/GSL-rl_COI_Folmer650pb_450taxa_1304seq.fasta")
DNA.folmer
```

    ## DNAStringSet object of length 1304:
    ##        width seq                                            names               
    ##    [1]   640 TGGATCGTTTGCTGCAATGGTA...ATTCTTTATCAACATTTATTT GBCI4815-14_Root_...
    ##    [2]   640 TGGATTGTTTGCTGCAATGGTA...ATTCTTTATCAACATTTATTT GBCI4819-14_Root_...
    ##    [3]   640 CGGATTGTATGCTGCAATGGTA...ATTCTTTATCAACATTTATTT GBCI4816-14_Root_...
    ##    [4]   640 TGGATTGTATGCTGCAATGGTA...ATTCTTTATCAACATTTATTT GBCI4820-14_Root_...
    ##    [5]   640 TGGATTGTTTGCTGCAATGGTA...ATTCTTTATCAACATTTATTT GBCI4817-14_Root_...
    ##    ...   ... ...
    ## [1300]   661 GACTCTTTATCTATATAGTGGG...GTCTTGTTTCAACATTTGTTC DESS022-20_Root_A...
    ## [1301]   661 GACTCTTTATCTATATAGTGGG...GTCTTGTTTCAACATTTGTTC STAN008-14_Root_A...
    ## [1302]   661 GACTCTTTATCTATATAGTGGG...GTCTTGTTTCAACATTTGTTC DESS018-20_Root_A...
    ## [1303]   658 AACCCTTTATCTGTATAGCGGA...CCCGTTTTGTTCCAACATTTG ECMOL120-11_Root_...
    ## [1304]   657 ACTCTTTACCTTTATAGAGGGG...CCTGTTCTGTTTCAACATTTA SERCI583-14_Root_...

Sequence name follows this nomenclature: *BOLD/NCBI Unique ID* \_ Root
\_ *Kingdom* \_ *Phylum* \_ *Class* \_ *Order* \_ *Family* \_ *Genus* \_
*Species*

``` r
names(DNA.folmer)[c(1, 10, 100)]
```

    ## [1] "GBCI4815-14_Root_Animalia_Cnidaria_Hydrozoa_Siphonophorae_Agalmatidae_Nanomia_Nanomia_cara"   
    ## [2] "BBPS066-19_Root_Animalia_Annelida_Polychaeta_Spionida_Spionidae_Laonice_Laonice_cirrata"      
    ## [3] "GBCM20402-19_Root_Animalia_Arthropoda_Hexanauplia_Calanoida_Temoridae_Temora_Temora_stylifera"

IDtaxa trained dataset

Metadata included the species rank category, that can be retrieve
easily.

``` r
metadata   <- readr::read_csv2("GSL-rl_COI/GSL-rl_COI_metadata.csv")

metadata %>% dplyr::filter(str_detect(Name, "Ammodytes")) %>% dplyr::select(Name, SEQavailable, BIN, SPsharingBIN, Validity) 
```

    ## # A tibble: 3 x 5
    ##   Name                 SEQavailable BIN          SPsharingBIN     Validity      
    ##   <chr>                <chr>        <chr>        <chr>            <chr>         
    ## 1 Ammodytes americanus Yes          BOLD:AAB4332 Ammodytes dubius Unreliable - ~
    ## 2 Ammodytes dubius     Yes          BOLD:AAB4332 Ammodytes ameri~ Unreliable - ~
    ## 3 Ammodytes hexapterus Yes          BOLD:AAB8000 <NA>             Reliable

``` r
metadata %>% dplyr::filter(family == "Asteriidae") %>% dplyr::select(Name, SEQavailable, BIN, SPsharingBIN, SPmissingSEQ.genus, Validity) 
```

    ## # A tibble: 9 x 6
    ##   Name       SEQavailable BIN        SPsharingBIN    SPmissingSEQ.ge~ Validity  
    ##   <chr>      <chr>        <chr>      <chr>           <chr>            <chr>     
    ## 1 Asterias ~ Yes          BOLD:AAB6~ <NA>            <NA>             Reliable  
    ## 2 Asterias ~ Yes          BOLD:AAB8~ <NA>            <NA>             Reliable  
    ## 3 Leptaster~ Yes          BOLD:AAB2~ Leptasterias l~ Leptasterias te~ Unreliabl~
    ## 4 Leptaster~ Yes          BOLD:AAB2~ Leptasterias g~ Leptasterias te~ Unreliabl~
    ## 5 Leptaster~ Yes          BOLD:AAJ1~ <NA>            Leptasterias te~ Unreliabl~
    ## 6 Leptaster~ Yes          BOLD:ACE6~ <NA>            Leptasterias te~ Unreliabl~
    ## 7 Leptaster~ No           <NA>       <NA>            <NA>             No sequen~
    ## 8 Stephanas~ Yes          BOLD:AAF3~ <NA>            <NA>             Reliable  
    ## 9 Urasteria~ Yes          BOLD:AAF5~ <NA>            <NA>             Reliable

The species rank category can be explore and add as a layer of
information on metabarcoding results.

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Requirements

*Optional section.* List the input data requirements or software
requirements to successfully execute the code.

## Caveats

Anything other users should be aware of including gaps in the input or
output data and warnings about appropriate use.

## Uncertainty

*Optional section.* Is there some uncertainty associated with the
output? Assumptions that were made?

## Acknowledgements

*Optional section.* List any contributors and acknowledge relevant
people or institutions

## References

*Optional section.*
