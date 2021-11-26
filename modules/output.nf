process OUTDIR {
	tag "outdir"

	input:
		path outdir

	output:
		path outdir

	when:
		params.init_outdir

	script:
	"""
	mkdir -p $init_outdir
	mkdir -p $init_outdir"/QC"
	mkdir -p $init_outdir"/REFERENCE"
	"""
}