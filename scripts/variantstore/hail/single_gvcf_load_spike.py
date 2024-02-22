import hail as hl
import os

TMP_DIR='gs://gvs-internal-scratch/mcovarr/scratch/hail_tmp'
OUTPUT_VDS_PREFIX='gs://gvs-internal-scratch/mcovarr/scratch/hail_combiner'

hl.init(tmp_dir=TMP_DIR, default_reference='GRCh38')

vcfs = (
    "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00405.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
    # "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00408.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
    # "gs://gvs-internal-quickstart/reblocked-v2-vcfs/HG00418.haplotypeCalls.er.raw.vcf.gz.rb.g.vcf.gz",
)

for vcf in vcfs:
    # strip the gs://bucket/prefixes... to get just the filename, then remove everything after the first dot
    id = os.path.basename(vcf).split('.')[0]
    combiner = hl.vds.new_combiner(
        output_path=f'{OUTPUT_VDS_PREFIX}/{id}.vds',
        temp_path=TMP_DIR,
        gvcf_paths=[vcf],
        use_genome_default_intervals=True,
    )
    combiner.run()
