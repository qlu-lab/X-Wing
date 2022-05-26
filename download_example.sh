#!/bin/bash

# X-Wing example data download #

## ---- Download the example data ---- ##
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/example/X-Wing_example.tar.gz
tar -zxvf X-Wing_example.tar.gz
rm -rf X-Wing_example.tar.gz

## ---- Download the required reference panel ---- ##
cd example/data
## ---- Step1: LOGODetect ---- ##
### EUR reference panel
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LOGODetect/LOGODetect_1kg_EUR.tar.gz
tar -zxvf LOGODetect_1kg_EUR.tar.gz
rm -rf LOGODetect_1kg_EUR.tar.gz

# EAS reference panel
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LOGODetect/LOGODetect_1kg_EAS.tar.gz
tar -zxvf LOGODetect_1kg_EAS.tar.gz
rm -rf LOGODetect_1kg_EAS.tar.gz

## ---- Step2: PANTHER ---- ##
### EUR pre-computed LD matrix
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/PANTHER/PANTHER_1kg_EUR.tar.gz
tar -zxvf PANTHER_1kg_EUR.tar.gz
rm -rf PANTHER_1kg_EUR.tar.gz

### EAS pre-computed LD matrix
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/PANTHER/PANTHER_1kg_EAS.tar.gz
tar -zxvf PANTHER_1kg_EAS.tar.gz
rm -rf PANTHER_1kg_EAS.tar.gz

### SNP information file
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/PANTHER/snpinfo_mult_1kg_hm3_PANTHER.tar.gz
tar -zxvf snpinfo_mult_1kg_hm3_PANTHER.tar.gz
rm -rf snpinfo_mult_1kg_hm3_PANTHER.tar.gz


## ---- Step3: LEOPARD ---- ##
### EUR reference panel
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LEOPARD/LEOPARD_1kg_hm3_EUR.tar.gz
tar -zxvf LEOPARD_1kg_hm3_EUR.tar.gz
rm -rf LEOPARD_1kg_hm3_EUR.tar.gz

### EUR pre-computed LD matrix
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LEOPARD/PANTHER_LEOPARD_1kg_EUR.tar.gz
tar -zxvf PANTHER_LEOPARD_1kg_EUR.tar.gz
rm -rf PANTHER_LEOPARD_1kg_EUR.tar.gz

### EAS reference panel
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LEOPARD/LEOPARD_1kg_hm3_EAS.tar.gz
tar -zxvf LEOPARD_1kg_hm3_EAS.tar.gz
rm -rf LEOPARD_1kg_hm3_EAS.tar.gz

### EAS pre-computed LD matrix
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LEOPARD/PANTHER_LEOPARD_1kg_EAS.tar.gz
tar -zxvf PANTHER_LEOPARD_1kg_EAS.tar.gz
rm -rf PANTHER_LEOPARD_1kg_EAS.tar.gz

### SNP information file
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LEOPARD/snpinfo_mult_1kg_hm3_PANTHER_LEOPARD.tar.gz
tar -zxvf snpinfo_mult_1kg_hm3_PANTHER_LEOPARD.tar.gz
rm -rf snpinfo_mult_1kg_hm3_PANTHER_LEOPARD.tar.gz
