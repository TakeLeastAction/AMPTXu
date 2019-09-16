#!/bin/bash
ifort -o coal CoalNuclAmptFinal.f90
cd submit_coal
condor_submit coal.con

watch condor_q lwchen
