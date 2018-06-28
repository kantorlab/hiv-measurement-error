import os

env = Environment(ENV=os.environ)
env.CacheDir("cache")

cpus = Value(8)

# Reference data

env.Command("scratch/pol.fa",
            ["lib/trim-reference.py",
             "data/HIV1_ALL_2016_2253-3870_DNA.fasta"],
            "python $SOURCES > $TARGET")

env.Command("scratch/pol.hmm",
            "scratch/pol.fa",
            "hmmbuild $TARGET $SOURCE")

# Downloads

env.Command("scratch/shiver-1.3.3.zip",
            Value("https://github.com/ChrisHIV/shiver/archive/v1.3.3.zip"),
            "wget -O $TARGET $SOURCE")

datasets = {"5VM":  "SRR961514",
            "PL11": "SRR6725661",
            "PL19": "SRR6725662"}

for dataset, accession in datasets.items():
    env.Command(["scratch/{}_1.fastq".format(accession),
                 "scratch/{}_2.fastq".format(accession)],
                Value(accession),
                " && ".join([
                     "cd scratch",
                     "prefetch $SOURCE",
                     "fastq-dump --defline-qual '+' --split-files --defline-seq '@$$sn[_$$rn]/$$ri' $SOURCE"]))

#datasets["PID"] = ["PID"]

# Software installs

#env.Command("scratch/shiver-install.log",
#            "lib/shiver-install.sh",
#            "$SOURCE > $TARGET")

# hivmmer

for dataset, accession in datasets.items():
    env.Command("scratch/{}.hmmsearch2.codons.csv".format(dataset),
                [Value("--id"), Value("scratch/{}".format(dataset)),
                 Value("--fq1"), "scratch/{}_1.fastq".format(accession),
                 Value("--fq2"), "scratch/{}_2.fastq".format(accession),
                 Value("--ref"), "scratch/pol.hmm",
                 Value("--cpu"), cpus],
                "hivmmer $SOURCES")

# hydra

#for dataset in datasets:
#    env.Command("scratch/hydra.{}".format(dataset),
#                ["scratch/{}_1.fastq".format(dataset),
#                 "scratch/{}_2.fastq".format(dataset)],
#                "quasitools hydra -o $TARGET $SOURCES")

# shiver

# bowtie2


# Results - codon tables

