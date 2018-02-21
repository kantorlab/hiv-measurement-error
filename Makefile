CPU = 16

# Reference

scratch/pol: lib/hmmbuild.sh data/ref.pol.fa
	$^ $@ &>$@.log
	touch $@

# Software installs

scratch/shiver: lib/shiver-install.sh
	$^ &>$@-install.log

scratch/PrimerID: lib/pid-install.sh
	$^ &>$@-install.log

# Data sets

scratch/5VM: lib/download.sh
	$< SRR961514 5VM &>$@.log
	touch $@

scratch/PL11: lib/download.sh
	$< SRR6725661 PL11 &>$@.log
	touch $@

scratch/PL19: lib/download.sh
	$< SRR6725662 PL19 &>$@.log
	touch $@

scratch/PID.download: lib/pid-download.sh
	$< &>$@.log
	touch $@

scratch/PID.link: lib/pid-link.sh scratch/PrimerID scratch/PID.download
	$< &>$@.log
	touch $@

scratch/PID.$(ID): lib/pid-align.sh scratch/PID.link
	$< $(ID) $(CPU) &>$@.log
	touch $@

scratch/PID: lib/pid-reads.sh scratch/PID.3223a scratch/PID.3223b scratch/PID.3223c scratch/PID.3236a scratch/PID.3236b scratch/PID.3236c
	$< &>$@.log
	touch $@

# Intermediate files

scratch/$(ID).iva: lib/iva.sh scratch/$(ID)
	$< $(ID) $(CPU) &>$@.log

scratch/$(ID).shiver: lib/shiver.sh scratch/shiver scratch/$(ID).iva scratch/pol
	$< $(ID) &>$@.log

scratch/$(ID).coverage: lib/coverage.sh data/hxb2.fa scratch/$(ID)
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

# Results

results/PID.hivmmer.codons.csv: lib/hivmmer.sh scratch/PID scratch/pol
	$^ $(CPU) scratch/PID.hmmsearch2.codons.csv &>$@.log
	mv scratch/PID.hmmsearch1.codons.csv $@

results/PID.hivmmer-ml.codons.csv: lib/model-codons.sh lib/model-codons.py scratch/model.PID scratch/model-cv
	$^ PID $@ &>$@.log

results/$(ID).hivmmer.codons.csv: lib/hivmmer.sh scratch/$(ID) scratch/pol
	$^ $(CPU) $@ &>$@.log

results/$(ID).shiver.codons.csv: lib/shiver-codons.sh lib/shiver-codons.py scratch/$(ID).shiver
	$< $(ID) $@ &>$@.log

results/$(ID).hydra.codons.csv: lib/hydra-codons.sh lib/hydra-codons.py scratch/hydra/$(ID)_1.align.bam scratch/hydra/$(ID)_2.align.bam
	$^ $@ &>$@.log

results/PID.pidalyse.codons.csv: lib/pid-codons.sh lib/pid-codons.py scratch/PID
	$^ $@ &>$@.log

results/PID.calls.csv: lib/calls.py data/PID.ref.fa results/PID.hivmmer.codons.csv results/PID.hydra.codons.csv results/PID.pidalyse.codons.csv
	python $< PID hivmmer hivmmer-ml hydra pidalyse &>$@.log

results/$(ID).calls.csv: lib/calls.py data/$(ID).ref.fa results/$(ID).hivmmer.codons.csv results/$(ID).hydra.codons.csv results/$(ID).shiver.codons.csv
	python $< $(ID) hivmmer hydra shiver &>$@.log

results/calls.stats.txt: lib/stat-calls.py results/5VM.calls.csv results/PL11.calls.csv results/PL19.calls.csv
	python $^ 1>$@ 2>$@.log

results/calls.eps: lib/plot-calls.py results/5VM.calls.csv results/PL11.calls.csv results/PL19.calls.csv
	python $< $@ &>$@.log

results/PID-detail.eps: lib/plot-detail.py data/PID.ref.fa results/PID.hivmmer.codons.csv results/PID.hydra.codons.csv results/PID.pidalyse.codons.csv
	python $< PID hivmmer hivmmer-ml hydra pidalyse &>$@.log

results/$(ID)-detail.eps: lib/plot-detail.py data/$(ID).ref.fa results/$(ID).hivmmer.codons.csv results/$(ID).hydra.codons.csv results/$(ID).shiver.codons.csv
	python $< $(ID) hivmmer hydra shiver &>$@.log

results/coverage.eps: lib/plot-coverage.py data/hxb2.fa scratch/5VM.coverage
	python $< &>$@.log

results/errors.eps: lib/plot-errors.py data/5VM.ref.fa data/PL11.ref.fa data/PL19.ref.fa results/5VM.hivmmer.codons.csv results/5VM.hydra.codons.csv results/5VM.shiver.codons.csv results/PL11.hivmmer.codons.csv results/PL11.hivmmer.codons.csv results/PL11.hydra.codons.csv results/PL11.shiver.codons.csv results/PL19.hivmmer.codons.csv results/PL19.hydra.codons.csv results/PL19.shiver.codons.csv
	python $< $@ &>$@.log

results/primer-id.eps: lib/plot-primer-id.py data/PID.ref.fa results/PID.hivmmer.codons.csv results/PID.hivmmer-ml.codons.csv results/PID.hydra.codons.csv results/PID.pidalyse.codons.csv
	python $< $@ &>$@.log

results/model.eps: lib/plot-model.py scratch/model-cv
	python $< $@ &>$@.log

