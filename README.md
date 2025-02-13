
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rtcassignment

<!-- badges: start -->

<!-- badges: end -->

The goal of rtcassignment is to assign output classes (Perfect ID, Genus
Group, New Taxonomic Group, etc) to microbial isolates based on
MALDI-TOF Real-Time Classification (RTC). This package takes the result
table in a text file format, applies predetermined logic for processing
those results and assignes classes, and tabulates the assignment.

What this program does:  
Step1: Accepts the raw data and cleans it up.  
Step2: Creates a flat file (each id on one line). This is the intermed
file object.  
Step3: Uses part of the intermediate file as final input, to generate
the result file, based on logic provided.  
Step4: Writes a text file as output with the results.

Computational Requirements Note this program is written for easy
implementation, not speed.  
Also it uses base R code, no libraries, so can be run on any machine
with R installed.

## Installation

You can install the development version of this package from
[GitHub](https://github.com/cm0109) with:

``` r
install.packages("devtools")
devtools::install_github("cm0109/rtcassignment")
```

## Example

This is a basic example which shows you how to use this package, using
the provided example data (“example\_data.txt”). The output will be
saved in a text file named “example\_assignment.txt” with this code.

``` r
library(rtcassignment)
rtc_to_results("inst/extdata/example_data.txt", "example_assignment.txt")
```
