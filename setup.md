# setup gatk
## download gatk zip file
https://github.com/broadinstitute/gatk/releases

## install java
https://adoptopenjdk.net/installation.html?variant=openjdk8&jvmVariant=hotspot#linux-pkg

## export path of gatk
export PATH="/path/to/gatk-package/:$PATH"

## set up gatk conda env
conda env create -f gatkcondaenv.yml
conda activate gatk


# setup vep

## clone vep from github
git clone https://github.com/Ensembl/ensembl-vep.git

## install dependencies
gcc
g++
make
perl -v

## install perl module
cpan App::cpanminus
cpanm Archive::Zip
cpanm DBD::mysql



cd ensembl-vep/


