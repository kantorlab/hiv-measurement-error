import os

env = Environment(ENV=os.environ,
                  LANG="en_US.utf-8",
                  LC_ALL="en_US.utf-8",
                  CPUS=1)

env.CacheDir("cache")
methods = ["hivmmer", "shiver", "hydra", "bowtie2", "bowtie2-pear"]
methods = {"5VM": methods, "PL11": methods, "PL19": methods, "PID": ["hivmmer", "pidalyse", "hydra", "bowtie2"]}

# Uncomment the second line to use the SLURM scheduler to run tasks in parallel.
# If your SLURM culster requires additional parameters, you can include them in
# this command.
srun = None
#srun = "srun"

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
            "data/hxb2-pol.fa",
            "bowtie2-build $SOURCE $TARGET > $TARGET")

# Downloads

env.Command("scratch/primer-id-5vm.zip",
            Value("https://ndownloader.figshare.com/articles/6713132/versions/1"),
            "wget -O $TARGET $SOURCE")

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

env.Command(["scratch/primer-id-5vm/3223a_R1.fastq.gz",
             "scratch/primer-id-5vm/3223a_R2.fastq.gz",
             "scratch/primer-id-5vm/3223b_R1.fastq.gz",
             "scratch/primer-id-5vm/3223b_R2.fastq.gz",
             "scratch/primer-id-5vm/3223c_R1.fastq.gz",
             "scratch/primer-id-5vm/3223c_R2.fastq.gz",
             "scratch/primer-id-5vm/3236a_R1.fastq.gz",
             "scratch/primer-id-5vm/3236a_R2.fastq.gz",
             "scratch/primer-id-5vm/3236b_R1.fastq.gz",
             "scratch/primer-id-5vm/3236b_R2.fastq.gz",
             "scratch/primer-id-5vm/3236c_R1.fastq.gz",
             "scratch/primer-id-5vm/3236c_R2.fastq.gz",
             "scratch/primer-id-5vm/3223a_QC_1_cons.fasta",
             "scratch/primer-id-5vm/3223b_QC_1_cons.fasta",
             "scratch/primer-id-5vm/3223c_QC_1_cons.fasta",
             "scratch/primer-id-5vm/3236a_QC_1_cons.fasta",
             "scratch/primer-id-5vm/3236b_QC_1_cons.fasta",
             "scratch/primer-id-5vm/3236c_QC_1_cons.fasta"],
            "scratch/primer-id-5vm.zip",
            "mkdir -p scratch/primer-id-5vm && unzip -o $SOURCE -d scratch/primer-id-5vm/")


for i in (1, 2):
    env.Command("scratch/PID_{}.fastq".format(i),
                ["scratch/primer-id-5vm/3223a_R{}.fastq.gz".format(i),
                 "scratch/primer-id-5vm/3223b_R{}.fastq.gz".format(i),
                 "scratch/primer-id-5vm/3223c_R{}.fastq.gz".format(i),
                 "scratch/primer-id-5vm/3236a_R{}.fastq.gz".format(i),
                 "scratch/primer-id-5vm/3236b_R{}.fastq.gz".format(i),
                 "scratch/primer-id-5vm/3236c_R{}.fastq.gz".format(i)],
                "zcat $SOURCES > $TARGET")

fastq["PID"] = ("scratch/PID_1.fastq", "scratch/PID_2.fastq")

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

    if dataset != "PID":

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
                rmdir=hydra_dir,
                mem_per_cpu=8)

    ## bowtie2

    SrunCommand("scratch/bowtie2.{}.bam".format(dataset),
                [Value("-x"), "scratch/pol.idx",
                 Value("-1"), fq1,
                 Value("-2"), fq2,
                 Value("--very-sensitive-local -X 1400")],
                "bowtie2 --threads $CPUS $SOURCES | samtools view -bS - > $TARGET",
                cpus=8)

    ## bowtie2-pear

    SrunCommand("scratch/bowtie2-pear.{}.bam".format(dataset),
                [Value("-x"), "scratch/pol.idx",
#                 Value("-1"), "scratch/{}.pear.unassembled.forward.fastq".format(dataset),
#                 Value("-2"), "scratch/{}.pear.unassembled.reverse.fastq".format(dataset),
                 Value("-U"), "scratch/{}.pear.assembled.fastq".format(dataset),
                 Value("--very-sensitive-local")],
                "bowtie2 --threads $CPUS $SOURCES | samtools view -bS - > $TARGET",
                cpus=8)

# Codon frequency tables

    ## hivmmer
    if dataset == "PID":
        # Use the first alignment for PID, since it is a sparse fragment with a
        # gap, and the second sample-specific alignment is no longer in HXB2
        # coordinates.
        hmmsearch = 1
    else:
        hmmsearch = 2

    env.Command("results/{}.hivmmer.codons.csv".format(dataset),
                "scratch/{}.hmmsearch{}.codons.csv".format(dataset, hmmsearch),
                "cp $SOURCE $TARGET")

    env.Command("results/{}.hivmmer.aavf".format(dataset),
                "scratch/{}.hmmsearch{}.aavf".format(dataset, hmmsearch),
                "cp $SOURCE $TARGET")

    ## shiver

    if dataset != "PID":
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
                 Value(0),
                 "scratch/bowtie2.{}.bam".format(dataset)],
                "python $SOURCES $TARGET")

    if dataset == "5VM":
        env.Command("results/{}.coverage.tsv".format(dataset),
                    "scratch/bowtie2.{}.bam".format(dataset),
                    "samtools view -f 3 $SOURCE | cut -f 4,9 > $TARGET")

    ## bowtie2-pear

    if dataset != "PID":
        SrunCommand("results/{}.bowtie2-pear.codons.csv".format(dataset),
                    ["lib/pileup.py",
                     Value(0),
                     "scratch/bowtie2-pear.{}.bam".format(dataset)],
                    "python $SOURCES $TARGET")

## pidalyse

env.Command("results/PID.pidalyse.codons.csv",
            ["lib/pid-codons.py",
             "scratch/primer-id-5vm/3223a_QC_1_cons.fasta",
             "scratch/primer-id-5vm/3223b_QC_1_cons.fasta",
             "scratch/primer-id-5vm/3223c_QC_1_cons.fasta",
             "scratch/primer-id-5vm/3236a_QC_1_cons.fasta",
             "scratch/primer-id-5vm/3236b_QC_1_cons.fasta",
             "scratch/primer-id-5vm/3236c_QC_1_cons.fasta"],
            "python $SOURCES > $TARGET")

# Variant distributions

for dataset in methods:

    env.Command("results/{}.variants.csv".format(dataset),
                ["lib/variants.py",
                 "data/{}.ref.fa".format(dataset),
                 Value(",".join(methods[dataset]))] + \
                ["results/{}.{}.codons.csv".format(dataset, method)
                 for method in methods[dataset]],
                "python $SOURCES $TARGET")

    env.Command("results/{}.variants.thresholded.csv".format(dataset),
                ["lib/variants-threshold.py",
                 "data/{}.ref.fa".format(dataset),
                 Value(",".join(methods[dataset]))] + \
                ["results/{}.{}.codons.csv".format(dataset, method)
                 for method in methods[dataset]],
                "python $SOURCES $TARGET")

    env.Command("results/{}.variants.test.log".format(dataset),
                ["lib/variants-test.py",
                 "results/{}.variants.thresholded.csv".format(dataset),
                 Value(",".join(methods[dataset]))],
                "python $SOURCES > $TARGET")

    env.Command("results/{}.top-errors.csv".format(dataset),
                ["lib/top-errors.py", "results/{}.variants.csv".format(dataset), Value(",".join(methods[dataset]))],
                "python $SOURCES $TARGET")

# Figures

    ## Supplementary 1-4

for dataset in methods:
    env.Command("results/{}-detail.eps".format(dataset),
                ["lib/plot-detail.py", Value(dataset), "data/{}.ref.fa".format(dataset)] + \
                ["results/{}.{}.codons.csv".format(dataset, method) for method in methods[dataset]],
                "python $SOURCES $TARGET")

# Figure 1

env.Command("results/coverage.eps",
            ["lib/plot-coverage.py", "results/5VM.coverage.tsv"],
            "python $SOURCES $TARGET")

# Figure 2

datasets = ("5VM", "PL11", "PL19")
env.Command("results/errors.eps",
            ["lib/plot-errors.py", Value(",".join(datasets)), Value(",".join(methods["5VM"]))] + \
            ["data/{}.ref.fa".format(dataset)
             for dataset in datasets] + \
            ["results/{}.{}.codons.csv".format(dataset, method)
             for method in methods["5VM"]
             for dataset in datasets],
            "python $SOURCES $TARGET")

# Figure 3

datasets = ("5VM", "PL11", "PL19")
env.Command("results/variant-dist.eps",
            ["lib/plot-variant-dist.py", Value(",".join(datasets)), Value(",".join(methods["5VM"]))] + \
            ["results/{}.variants.thresholded.csv".format(dataset) for dataset in datasets],
            "python $SOURCES $TARGET")

# Figure 4

env.Command("results/primer-id.eps",
            ["lib/plot-primer-id.py", "data/PID.ref.fa", "results/PID.variants.thresholded.csv"] + \
            ["results/PID.{}.codons.csv".format(method) for method in methods["PID"]],
            "python $SOURCES $TARGET")

# Read counts

for dataset in methods:
    SrunCommand("results/{}.hivmmer_read_count.log".format(dataset),
                ["lib/hivmmer-count.py",
                 "scratch/{}.hmmsearch2.txt".format(dataset)],
                "python $SOURCES > $TARGET",
                mem_per_cpu=8)
    SrunCommand("results/{}.hydra_read_count.log".format(dataset),
                "scratch/hydra.{}/align.bam".format(dataset),
                "samtools flagstat $SOURCE > $TARGET")
    SrunCommand("results/{}.bowtie2_read_count.log".format(dataset),
                "scratch/bowtie2.{}.bam".format(dataset),
                "samtools flagstat $SOURCE > $TARGET")

for dataset in datasets:
    SrunCommand("results/{}.shiver_read_count.log".format(dataset),
                "scratch/shiver.{0}/{0}_remap.bam".format(dataset),
                "samtools flagstat $SOURCE > $TARGET")
    SrunCommand("results/{}.bowtie2-pear_read_count.log".format(dataset),
                "scratch/bowtie2-pear.{}.bam".format(dataset),
                "samtools flagstat $SOURCE > $TARGET")

# vim: syntax=python expandtab sw=4 ts=4
