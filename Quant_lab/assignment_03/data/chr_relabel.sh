#!/bin/bash

awk '{gsub(/chrXVI/, "16"); 
      gsub(/chrXV/, "15"); 
      gsub(/chrXIV/, "14");
      gsub(/chrXIII/, "13"); 
      gsub(/chrXII/, "12"); 
      gsub(/chrXI/, "11"); 
      gsub(/chrX/, "10"); 
      gsub(/chrIX/, "9"); 
      gsub(/chrVIII/, "8");
      gsub(/chrVII/, "7"); 
      gsub(/chrVI/, "6"); 
      gsub(/chrV/, "5");
      gsub(/chrIV/, "4"); 
      gsub(/chrIII/, "3"); 
      gsub(/chrII/, "2");
      gsub(/chrI/, "1");
      gsub(/chrM/, "M");
      print;}' /Users/cmdb/qbb2020-answers/Quant_lab/assignment_03/data/BYxRM_segs_saccer3.bam.simplified.vcf.copy > renum_chr.vcf