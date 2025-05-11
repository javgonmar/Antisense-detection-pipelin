# Antisense-detection-pipelin
The workflow design is primarily intended for the detection and quantification of antisense transcripts from any organism using raw ssRNA-seq data (for any number of conditions and replicates), based on Bash and R as programming engines. 

User Manual for the Pipeline

To define parameters and file paths, it is necessary to create a parameters file (params file). A template parameter file is provided to the user for completion. For each variable added, the structure shown before the colon (:) must be strictly followed. If additional variables not present in the template need to be added, their format must match that of the referenced examples, corresponding to the same type of file. The only variation allowed is in the first number (referring to the condition) and the second number (referring to the replicate). The numerical order of samples, replicates, and conditions must be respected. There is no limit to the number of conditions or samples that can be added.

The variables included in the parameter file are as follows (as they appear in the template):

number_conditions: Refers to the total number (in numeric format) of experimental conditions in the study. The expected input is a discrete numeric value.

number_fasta_condition: Refers to the number of FASTA files provided for each condition. The expected input is a discrete numeric value.

fasta_path_condition_<>_<>: Refers to the location of each FASTA file. The first number in the variable name refers to the condition number, and the second to the replicate number within that condition. The expected input is the absolute path to the desired FASTA file.

genome_path: Refers to the path of the genome file of the organism under study. The expected input is the absolute path to the genome FASTA file.

annotation_path: Refers to the path of the reference annotation file of the organism. The expected input is the absolute path to the GTF annotation file.

working_directory: Refers to the location where the project directory will be generated. The expected input is the absolute path to the desired terminal environment location for storing the output files.

installation_directory: Refers to the location where the pipeline scripts are installed on the system. The expected input is the absolute path to the script directory.

experiment_name: Desired name for the generated working directory.

number_samples: Refers to the total number of FASTA files provided, considering only the sample files (excluding the genome FASTA file). This must equal the sum of FASTA files across all conditions.

genome_parameter: Refers to the value associated with the genome size of the organism under study. The expected input is a discrete number or, at most, a number with one decimal place. It must be calculated using the following expression (where L is the genome length in base pairs):

\text{genome_parameter} = \min(14,\ \log_2((L / 2) - 1))
Once the parameter file has been properly completed, the pipeline can be executed from a terminal (with Bash access) using the following command:

bash
Copiar
Editar
bash antipipe.sh params.txt
The results will be organized into a new directory named after the specified experiment name. This directory will contain five subfolders: “genome”, “annotation”, “samples”, “scripts”, and “results”.

The "genome" and "annotation" folders will store files related to the genome and its processing.

The "scripts" folder will contain auto-generated elements required for the workflow to function properly.

The "samples" folder will organize the different sample files by condition, including the transcript quantification files for each replicate.

Finally, the "results" folder will include the new antisense annotation file along with other files related to its processing.
