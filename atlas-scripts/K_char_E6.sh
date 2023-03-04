#!/bin/bash


{ time ./K_char.py -g E6_s -f facetsE6dim0 -o e6dim0 -c 10 >e6_log;} 2> e6_log
{ time ./K_char.py -g E6_s -f facetsE6dim1 -o e6dim1 -c 10 >>e6_log;} 2>> e6_log
{ time ./K_char.py -g E6_s -f facetsE6dim2 -o e6dim2 -c 10 >>e6_log;} 2>> e6_log
{ time ./K_char.py -g E6_s -f facetsE6dim3 -o e6dim3 -c 10 >>e6_log;} 2>> e6_log
{ time ./K_char.py -g E6_s -f facetsE6dim4 -o e6dim4 -c 10 >>e6_log;} 2>> e6_log
{ time ./K_char.py -g E6_s -f facetsE6dim5 -o e6dim5 -c 10 >>e6_log;} 2>> e6_log
{ time ./K_char.py -g E6_s -f facetsE6dim6 -o e6dim6 -c 10 >>e6_log;} 2>> e6_log
wc -l E6_s/e6dim0* >> e6_log
wc -l E6_s/e6dim1* >> e6_log
wc -l E6_s/e6dim2* >> e6_log
wc -l E6_s/e6dim3* >> e6_log
wc -l E6_s/e6dim4* >> e6_log
wc -l E6_s/e6dim5* >> e6_log
wc -l E6_s/e6dim6* >> e6_log




