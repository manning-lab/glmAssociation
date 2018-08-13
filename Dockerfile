FROM robbyjo/r-mkl-bioconductor:3.4.3-16.04-2018.1

MAINTAINER tmajaria@broadinstitute.org

RUN apt-get update && apt-get -y install git dstat

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('data.table')"

RUN git clone https://github.com/manning-lab/glmAssociation.git && cd ./glmAssociation && git pull origin master