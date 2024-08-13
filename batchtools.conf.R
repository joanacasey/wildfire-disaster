default.resources = list(
  queue="mesa.q", 
  group="mesa", 
  ncpus=24, 
  blas.threads=1
)
cluster.functions = makeClusterFunctionsSGE(file.path(getwd(),"sge.tmpl")) 
