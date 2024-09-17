# Copyright (C) 2024 Troels Mouritzen
# 
# Author: Troels Mouritzen
# Affiliation: Aarhus University
# Email: twm@mbg.au.dk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# Description: This script takes a window-based file and performs logistic regression
# for to predict a binary outcome based on a null model and a model of interest.

# Associated publication: [PAPER TITLE], [JOURNAL], [YEAR]
# DOI: [DOI OF PAPER]

################################################################################
# Load packages and files
################################################################################

library(tidyverse)
windowfile <- read_csv("20240618_Suppl_File_1_windows.csv.gz")

################################################################################
# Logistic regression
################################################################################

dependent_variable <- "Bacterial_GWAS" # Which variable to predict
independent_variables_of_interest <- c("HAPk") # Which variables to include in the model
independent_variables_control <- c("GC") # Which variables to control for


var_int <- str_c(independent_variables_of_interest, collapse = " + ")
var_contr <- str_c(independent_variables_control, collapse = " + ")

your_model = as.formula(str_c(dependent_variable, 
                              " ~ ", var_int, " + ", var_contr))
your_null_model = as.formula(str_c(dependent_variable, " ~ ",
                                   var_contr))

glm(your_model, data = windowfile, family = "binomial") -> model
glm(your_null_model, data = windowfile, family = "binomial") -> null

anova(model, null, test = "Chisq")
summary(model)
