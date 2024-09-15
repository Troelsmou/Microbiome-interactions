# Mantel-test

This is the code for doing the mantel test for the paper.

## Software Requirements

The needed packages are listed in requirements.txt

## Installation

```bash
git clone (https://github.com/Troelsmou/Microbiome-interactions/edit/main/Mantel-test)

```

## Usage

The three scripts are supposed to be run interactively in an R session.

## Input files

`input_file.txt` [Description]

## Output files

`output_file.txt` [Description]

## Intermediate files

`intermediate_file.tmp` [Description]

## File relationships

@startuml
file "Input File 1" as input1
file "Input File 2" as input2
file "Shared File" as shared

rectangle "Layer 1" {
  [Process A] as A
  [Process B] as B
}

rectangle "Layer 2" {
  [Process C] as C
  [Process D] as D
}

input1 --> A
input2 --> B
shared --> A
shared --> B
shared --> C
A --> C
B --> D
@enduml

## Runtime
Typical runtime: X minutes
