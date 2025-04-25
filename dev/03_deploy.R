devtools::check()
devtools::build(path = getwd(), binary = T)

# add docker
golem::add_dockerfile_with_renv(output_dir = "docker", from = "rocker/shiny-verse")

# then modify the `Dockerfile` and `Dockerfile_base` to
# 1. add `./docker/` to the `COPY` command
# 2. modify the `Dockerfile`: FROM jinlongru/xxx_base:latest
