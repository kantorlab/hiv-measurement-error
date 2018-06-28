CPU = 8

# Reference data

scratch/pol: lib/hmmbuild.sh data/ref.pol.fa
	$^ $@ &>$@.log
	touch $@

# Software installs

scratch/shiver: lib/shiver-install.sh
	$^ &>$@-install.log

scratch/PrimerID: lib/pid-install.sh
	$^ &>$@-install.log

# Data sets

scratch/data.5VM: lib/download.sh
	$< SRR961514 5VM &>$@.log
	touch $@

scratch/data.PL11: lib/download.sh
	$< SRR6725661 PL11 &>$@.log
	touch $@

scratch/data.PL19: lib/download.sh
	$< SRR6725662 PL19 &>$@.log
	touch $@

scratch/download.PID: lib/pid-download.sh
	$< &>$@.log
	touch $@

scratch/link.PID: lib/pid-link.sh scratch/PrimerID scratch/download.PID
	$< &>$@.log
	touch $@

scratch/pidalyse.$(ID): lib/pid-align.sh scratch/link.PID scratch/PrimerID
	$< $(ID) $(CPU) &>$@.log
	touch $@

scratch/data.PID: lib/pid-reads.sh scratch/pidalyse.3223a scratch/pidalyse.3223b scratch/pidalyse.3223c scratch/pidalyse.3236a scratch/pidalyse.3236b scratch/pidalyse.3236c
	$< &>$@.log
	touch $@

# Intermediate files

scratch/hydra.$(ID): lib/hydra.sh scratch/data.$(ID)
	$< $(ID) $(CPU) &>$@.log

scratch/iva.$(ID): lib/iva.sh scratch/data.$(ID)
	$< $(ID) $(CPU) &>$@.log

scratch/shiver.$(ID): lib/shiver.sh scratch/shiver scratch/iva.$(ID) scratch/pol
	$< $(ID) &>$@.log

scratch/bowtie2.$(ID): lib/bowtie2.sh lib/consensus.py data/hxb2.fa scratch/data.${ID}
	$< $(ID) data/hxb2.fa $(CPU) &>$@.log
	touch $@

scratch/bowtie2-pear.$(ID): lib/bowtie2-pear.sh lib/consensus.py data/hxb2.fa scratch/data.${ID}
	$< $(ID) data/hxb2.fa $(CPU) &>$@.log
	touch $@

scratch/coverage.$(ID): lib/coverage.sh data/hxb2.fa scratch/data.$(ID)
	$< $(ID) data/hxb2.fa $(CPU) &>$@.log
	touch $@

scratch/model.PID: lib/model-prep.sh lib/model-prep.py results/PID.hivmmer.codons.csv
	$< lib/model-prep.py hmmsearch1 PID &>$@.log
	touch $@

scratch/model.$(ID): lib/model-prep.sh lib/model-prep.py results/$(ID).hivmmer.codons.csv
	$< lib/model-prep.py hmmsearch2 $(ID) &>$@.log
	touch $@

scratch/model-cv: lib/model-cv.py scratch/model.5VM scratch/model.PL11 scratch/model.PL19
	python $< $(CPU) &>$@.log
	touch $@

# Results - codon tables

results/PID.hivmmer.codons.csv: lib/hivmmer.sh scratch/data.PID scratch/pol
	$< PID scratch/pol $(CPU) scratch/PID.hmmsearch2.codons.csv &>$@.log
	mv scratch/PID.hmmsearch1.codons.csv $@

results/PID.hivmmer-ml.codons.csv: lib/model-codons.sh lib/model-codons.py scratch/model.PID scratch/model-cv
	$^ PID $@ &>$@.log

results/$(ID).hivmmer.codons.csv: lib/hivmmer.sh scratch/data.$(ID) scratch/pol
	$< $(ID) scratch/pol $(CPU) $@ &>$@.log

results/$(ID).shiver.codons.csv: lib/shiver-codons.sh lib/pileup.py scratch/shiver.$(ID)
	$< $(ID) $@ &>$@.log

results/$(ID).hydra.codons.csv: lib/hydra-codons.sh lib/hydra-codons.py scratch/hydra.$(ID)
	$^ $@ &>$@.log

results/$(ID).bowtie2.codons.csv: lib/bowtie2-codons.sh lib/pileup.py data/ref.pol.fa scratch/bowtie2.$(ID)
	$< $(ID) data/ref.pol.fa $@ &>$@.log

results/$(ID).bowtie2-pear.codons.csv: lib/bowtie2-pear-codons.sh lib/pileup.py data/ref.pol.fa scratch/bowtie2-pear.$(ID)
	$< $(ID) data/ref.pol.fa $@ &>$@.log

results/PID.pidalyse.codons.csv: lib/pid-codons.sh lib/pid-codons.py scratch/data.PID
	$^ $@ &>$@.log

# Results - call distributions

results/PID.calls.csv: lib/calls.py data/PID.ref.fa results/PID.hivmmer.codons.csv results/PID.hydra.codons.csv results/PID.pidalyse.codons.csv
	python $< PID hivmmer hivmmer-ml hydra pidalyse &>$@.log

results/$(ID).calls.csv: lib/calls.py data/$(ID).ref.fa results/$(ID).hivmmer.codons.csv results/$(ID).hydra.codons.csv results/$(ID).shiver.codons.csv
	python $< $(ID) hivmmer hydra shiver &>$@.log

results/calls.stats.txt: lib/stat-calls.py results/5VM.calls.csv results/PL11.calls.csv results/PL19.calls.csv
	python $^ 1>$@ 2>$@.log

results/calls.eps: lib/plot-calls.py results/5VM.calls.csv results/PL11.calls.csv results/PL19.calls.csv
	python $< $@ &>$@.log

# Results - alignment details

results/PID-detail.eps: lib/plot-detail.py data/PID.ref.fa results/PID.hivmmer.codons.csv results/PID.hydra.codons.csv results/PID.pidalyse.codons.csv
	python $< PID hivmmer hivmmer-ml hydra pidalyse &>$@.log

results/$(ID)-detail.eps: lib/plot-detail.py data/$(ID).ref.fa results/$(ID).hivmmer.codons.csv results/$(ID).hydra.codons.csv results/$(ID).shiver.codons.csv results/$(ID).bowtie2.codons.csv results/$(ID).bowtie2-pear.codons.csv
	python $< $(ID) hivmmer hydra shiver bowtie2 bowtie2-pear &>$@.log

# Results - remaining plots

results/coverage.eps: lib/plot-coverage.py data/hxb2.fa scratch/coverage.5VM
	python $< &>$@.log

results/errors.eps: lib/plot-errors.py data/5VM.ref.fa data/PL11.ref.fa data/PL19.ref.fa results/5VM.hivmmer.codons.csv results/5VM.hydra.codons.csv results/5VM.shiver.codons.csv results/5VM.bowtie2.codons.csv results/5VM.bowtie2-pear.codons.csv results/PL11.hivmmer.codons.csv results/PL11.hivmmer.codons.csv results/PL11.hydra.codons.csv results/PL11.shiver.codons.csv results/PL11.bowtie2.codons.csv results/PL11.bowtie2-pear.codons.csv results/PL19.hivmmer.codons.csv results/PL19.hydra.codons.csv results/PL19.shiver.codons.csv results/PL19.bowtie2.codons.csv results/PL19.bowtie2-pear.codons.csv
	python $< $@ &>$@.log

results/primer-id.eps: lib/plot-primer-id.py data/PID.ref.fa results/PID.hivmmer.codons.csv results/PID.hivmmer-ml.codons.csv results/PID.hydra.codons.csv results/PID.pidalyse.codons.csv
	python $< $@ &>$@.log

results/model.eps: lib/plot-model.py scratch/model-cv
	python $< $@ &>$@.log

