clc;
clear;
fclose all;

CloneCNA('./example_data/tumor.count','./example_data/normal.count', './example_data/tumor.depth', './example_data/target.gc', './results/', './plots/');