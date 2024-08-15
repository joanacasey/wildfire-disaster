Raw population data can be downloaded using `./code/get_data.sh`. 

On Brain, use the following modules to run code. 

```
module load 'R/R-4.3.0_GCC-10.3.0'
module load UDUNITS/udunits-2.2.26
module load GDAL/gdal-3.0.2 
module load GEOS/geos-3.8.0
module load PROJ/proj-6.2.1
module load GCC/gcc-10.3.0
module load pandoc/pandoc-2.7.2
```

Then use `Rscript code/kde_pop_criteria.r`. 

_Currently, preprocessing the full 30m CONUS population runs forever or runs out of memory and crashes._ 

* Minimally modified file from Joan is at `code/UTM_Zone_All_Pop_Density_Batch.r`. 

* Best attempt so far at `code/kde_pop_criteria.r`. 

* `batchtools.conf.R` and `sge.tmp` used to configure SGE on Brain cluster. 