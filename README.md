Setup:

1. Edit 'Snakefile' to include the appropriate rules. Most likely this is already set up correctly. 'get_final_targets' generates alignment files for all samples in the file 'targets.txt'. 'get_final_targets_all' generates alignment files for all samples included in the metadata tsv. 'get_multisample_only' gets only multisample vcfs and will need the commented lines uncommented to work.

2. Edit 'config/project_config.tsv' to include your sample data.
3. Edit 'config/targets.txt' to include only the samples you wish to generate alignments for.
4. Update 'config/config.yaml'
    - threads here is per sample.
5. Copy this directory structure into the final directory on /n and move into that directory
6. Run 'snakemake --use-conda --cores N' where N=the number of cores you would like to devote to the project. You should run this command within a minimal snakemake conda environment. The only other tool that snakemake expects to be available from your PATH is samtools (>=1.12 is fine). This will be changed in the future to require nothing but the minimal snakemake environment.