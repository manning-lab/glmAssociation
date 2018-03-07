FROM uwgac/r343-topmed:master
RUN apt-get update
RUN apt-get -y install git

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('data.table')"

RUN git clone https://github.com/manning-lab/glmAssociation.git && cd ./glmAssociation && git pull origin master