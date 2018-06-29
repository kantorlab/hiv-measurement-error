import os

env = Environment(ENV=os.environ,
                  LANG="en_US.utf-8",
                  LC_ALL="en_US.utf-8",
                  CPUS=1)

env.CacheDir("cache")

# Use the SLURM scheduler to run jobs by default. If you do not have SLURM, set
# the command to None to disable.
srun = "srun"
#srun = None

def SrunCommand(targets, sources, cmd, rmdir="", cpus=1, mem_per_cpu=3, timelimit="24:00:00"):
    global srun, env
    if srun is not None:
        cmd = "{} -c {} --mem={}G -t {} {}".format(srun, cpus, cpus*mem_per_cpu,
                                                   timelimit, cmd.replace("$CPUS", str(cpus)))
    if rmdir:
        assert not rmdir.startswith("/")
        cmd = "rm -rf {} && {}".format(rmdir, cmd)
    return env.Command(targets, sources, cmd)

# Reference data

env.Command("scratch/pol.fa",
            ["lib/trim-reference.py",
             "data/HIV1_ALL_2016_2253-3870_DNA.fasta"],
            "python $SOURCES > $TARGET")

env.Command("scratch/pol.hmm",
            "scratch/pol.fa",
            "hmmbuild $TARGET $SOURCE")

env.Command("scratch/pol.idx",
            "data/hxb2.fa",
            "bowtie2-build $SOURCE $TARGET > $TARGET")

# Downloads

env.Command("scratch/shiver-1.4.1.zip",
            Value("https://github.com/ChrisHIV/shiver/archive/v1.4.1.zip"),
            "wget -O $TARGET $SOURCE")

accessions = {"5VM": "SRR961514",
              "PL11": "SRR6725661",
              "PL19": "SRR6725662"}

for dataset, accession in accessions.items():
    env.Command(["scratch/{}_1.fastq".format(accession),
                 "scratch/{}_2.fastq".format(accession)],
                Value(accession),
                " && ".join([
                     "cd scratch",
                     "prefetch $SOURCE",
                     "fastq-dump --defline-qual '+' --split-files --defline-seq '@$$sn[_$$rn]/$$ri' $SOURCE"]))

fastq = dict((dataset, ("scratch/{}_1.fastq".format(accession), "scratch/{}_2.fastq".format(accession)))
             for dataset, accession in accessions.items())

#datasets["PID"] = ["PID"]

# Software installs

env.Command(["scratch/shiver/shiver_init.sh",
             "scratch/shiver/shiver_align_contigs.sh",
             "scratch/shiver/shiver_map_reads.sh",
             "scratch/shiver/config.sh"],
            ["lib/shiver-install.sh",
             "scratch/shiver-1.4.1.zip"],
            "bash $SOURCES")

SrunCommand("scratch/shiver/reference.log",
            ["scratch/shiver/shiver_init.sh",
             Value("scratch/shiver/reference"),
             "scratch/shiver/config.sh",
             "data/HIV1_ALL_2016_2253-3870_DNA.fasta",
             "data/adapters.fa",
             "data/primers.fa"],
            "bash $SOURCES > $TARGET")

# Alignments

for dataset, (fq1, fq2) in fastq.items():

    ## hivmmer

    SrunCommand(["scratch/{}.hmmsearch1.codons.csv".format(dataset),
                 "scratch/{}.hmmsearch1.aavf".format(dataset),
                 "scratch/{}.hmmsearch2.codons.csv".format(dataset),
                 "scratch/{}.hmmsearch2.aavf".format(dataset),
                 "scratch/{}.pear.assembled.fastq".format(dataset),
                 "scratch/{}.pear.unassembled.forward.fastq".format(dataset),
                 "scratch/{}.pear.unassembled.reverse.fastq".format(dataset),
                 "scratch/{}.pear.discarded.fastq".format(dataset)],
                [Value("--id"), Value("scratch/{}".format(dataset)),
                 Value("--fq1"), fq1,
                 Value("--fq2"), fq2,
                 Value("--ref"), "scratch/pol.hmm"],
                "hivmmer --cpu $CPUS $SOURCES",
                cpus=4)

    ## shiver

    iva_dir = "scratch/iva.{}".format(dataset)

    SrunCommand("{}/contigs.fasta".format(iva_dir),
                [Value("-f"), fq1, Value("-r"), fq2],
                "iva -t $CPUS $SOURCES {}".format(iva_dir),
                rmdir=iva_dir,
                cpus=8)

    shiver_dir = "scratch/shiver.{}".format(dataset)

    SrunCommand(["{}/{}_cut_wRefs.fasta".format(shiver_dir, dataset),
                 "{}/{}_raw_wRefs.fasta".format(shiver_dir, dataset),
                 "{}/{}.blast".format(shiver_dir, dataset)],
                ["lib/shiver-contigs.sh",
                 Value(shiver_dir),
                 "scratch/shiver/shiver_align_contigs.sh",
                 "scratch/shiver/reference.log",
                 "scratch/shiver/config.sh",
                 "{}/contigs.fasta".format(iva_dir),
                 Value(dataset)],
                "bash $SOURCES",
                rmdir=shiver_dir)

    SrunCommand(["{}/{}_remap.bam".format(shiver_dir, dataset),
                 "{}/{}_remap_ref.fasta".format(shiver_dir, dataset)],
                ["lib/shiver-reads.sh",
                 Value(shiver_dir),
                 "scratch/shiver/shiver_map_reads.sh",
                 "scratch/shiver/reference.log",
                 "scratch/shiver/config.sh",
                 "{}/contigs.fasta".format(iva_dir),
                 Value(dataset),
                 "{}/{}.blast".format(shiver_dir, dataset),
                 "{}/{}_cut_wRefs.fasta".format(shiver_dir, dataset),
                 fq1,
                 fq2],
                "bash $SOURCES")

    ## hydra

    hydra_dir = "scratch/hydra.{}".format(dataset)

    SrunCommand("{}/align.bam".format(hydra_dir),
                [fq1, fq2],
                "quasitools hydra -o {} $SOURCES".format(hydra_dir),
                rmdir=hydra_dir)

    ## bowtie2

    SrunCommand("scratch/bowtie2.{}.bam".format(dataset),
                [Value("-x"), "scratch/pol.idx",
                 Value("-1"), fq1,
                 Value("-2"), fq2],
                "bowtie2 --no-mixed --no-discordant --threads $CPUS $SOURCES | samtools view -bS - > $TARGET",
                cpus=8)

    ## bowtie2-pear

    SrunCommand("scratch/bowtie2-pear.{}.bam".format(dataset),
                [Value("-x"), "scratch/pol.idx",
                 Value("-1"), "scratch/{}.pear.unassembled.forward.fastq".format(dataset),
                 Value("-2"), "scratch/{}.pear.unassembled.reverse.fastq".format(dataset),
                 Value("-U"), "scratch/{}.pear.assembled.fastq".format(dataset)],
                "bowtie2 --no-mixed --no-discordant --threads $CPUS $SOURCES | samtools view -bS - > $TARGET",
                cpus=8)

# Results - codon tables

    ## hivmmer

    env.Command("results/{}.hivmmer.codons.csv".format(dataset),
                "scratch/{}.hmmsearch2.codons.csv".format(dataset),
                "cp $SOURCE $TARGET")

    env.Command("results/{}.hivmmer.aavf".format(dataset),
                "scratch/{}.hmmsearch2.aavf".format(dataset),
                "cp $SOURCE $TARGET")

    ## shiver

    SrunCommand("results/{}.shiver.codons.csv".format(dataset),
                ["lib/pileup.py",
                 Value(0),
                 "{}/{}_remap.bam".format(shiver_dir, dataset)],
                "python $SOURCES $TARGET")

    ## hydra

    SrunCommand("results/{}.hydra.codons.csv".format(dataset),
                ["lib/pileup.py",
                 Value(0),
                 "{}/align.bam".format(hydra_dir)],
                "python $SOURCES $TARGET")

    ## bowtie2

    SrunCommand("results/{}.bowtie2.codons.csv".format(dataset),
                ["lib/pileup.py",
                 Value(2252),
                 "scratch/bowtie2.{}.bam".format(dataset)],
                "python $SOURCES $TARGET")

    ## bowtie2-pear

    SrunCommand("results/{}.bowtie2-pear.codons.csv".format(dataset),
                ["lib/pileup.py",
                 Value(2252),
                 "scratch/bowtie2-pear.{}.bam".format(dataset)],
                "python $SOURCES $TARGET")

# vim: syntax=python expandtab sw=4 ts=4
