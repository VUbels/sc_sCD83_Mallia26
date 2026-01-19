setup_py_env <- function(py_env_name, py_location) {

options(reticulate.conda_binary = py_location) 
  
# Check existing environments
reticulate::conda_list()

# Create if doesn't exist, otherwise just use it
if (!(py_env_name %in% reticulate::conda_list()$name)) {
  # Environment doesn't exist - create it
  conda_create(py_env_name,
               python_version = "3.12",
               packages = c("pip", "umap-learn"))
  use_condaenv(paste0("/home/uvictor/miniconda3/envs/", py_env_name), required = TRUE)
  
  reticulate::py_install(
    packages = c("torch", "torchvision"),
    pip = TRUE,
    extra_options = c(
      "--index-url", "https://download.pytorch.org/whl/rocm6.4"
    )
  )
  
  reticulate::py_install(
    packages = "git+https://github.com/broadinstitute/CellBender.git@refs/pull/420/head",
    pip = TRUE
  )
  
} else {
  # Environment exists
  use_condaenv(paste0("/home/uvictor/miniconda3/envs/", py_env_name), required = TRUE)
}

# Verify configuration
py_config()

}