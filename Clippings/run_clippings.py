"""
Script to run both elements of clippings from a single command
"""

import argparse
import sys
import os
import logging
import glob
import datetime

from Clippings.count_miRNAs import main as count_miRNAs
from Clippings.deg_count_with_UMIs import main as deg_count_with_UMIs

## copied from deg_count_with_UMIs
def _parse_cmdl(cmdl):
    """ Define and parse command-line interface. """

    parser = argparse.ArgumentParser(
        description="Degradation and/or miRNA counting from TSO containing reads in a Cell Ranger outs folder "
                    "implementation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "CRouts", help="Path to Cellranger outs folder.")

    parser.add_argument(
        "-C", "--cores", default=10, help="Number of cores.")

    parser.add_argument(
        "-TSS", "--TSSgtf", help="Path to gtf file.", required=False, default=None)

    parser.add_argument('--outdir', dest='outdir', help='Output folder', default="Clippings")
    
    parser.add_argument('--genome', dest='genome',
                        help='Genome version to record in h5 file. eg. \'hg38\' or \'mm10\'', default=None)
   
    parser.add_argument('-miRNAgff3','--miRNAgff3', help='miRBase gff3 file', required=False, default=None)
    parser.add_argument('--write_bam_output',dest='write_bam_output',required=False, default=True,
                        help = 'Write candidate miRNA processing product reads/degraded reads to BAM output')

    return parser.parse_args(cmdl)

def main():
    """
    Break CLI arguments down into the  args needed by deg_count and count_miRNAs and launch them
    """
    cmdl=sys.argv[1:]
    args = _parse_cmdl(cmdl)

    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S",
                    level=logging.INFO, force=True)
    
    logging.info(f"Args: {args}")

# Handle different types of cellranger outs folders and check to make sure they are complete

    cellranger_folder_bam = str(args.CRouts)+"/possorted_genome_bam.bam"
    cellranger_folder_raw = str(args.CRouts)+"/raw_feature_bc_matrix/"

    not_10x_folder_error = """Not a proper 10x outs folder. 
                      A 10x outs folder must include one bam file labelled possorted_genome_bam.bam or
                        a contain a single bam file of the form sample_name_possorted_genome_bam.bam
                        with all other files having the same prefix and a raw matrix directory. 
                        If the raw matrix is tarballed unpack before running Clippings"""

    cellranger_folder_bam = str(args.CRouts)+"/possorted_genome_bam.bam"
    cellranger_multi_folder_bam = str(args.CRouts)+"/sample_alignments.bam"
    cellranger_folder_raw = str(args.CRouts)+"/raw_feature_bc_matrix/"
    cellranger_multi_folder_filtered = str(args.CRouts)+"/sample_feature_bc_matrix/"

    # Check if the "possorted_genome_bam.bam" file exists
    if not os.path.isfile(cellranger_folder_bam):
        # Check if there is only one file in the directory with the pattern "*possorted_genome_bam.bam"
        bam_files = glob.glob(str(args.CRouts) + "/*possorted_genome_bam.bam")
        if len(bam_files) == 1:
            cellranger_folder_bam = bam_files[0]
            cellranger_folder_raw = bam_files[0].replace("possorted_genome_bam.bam", "raw_feature_bc_matrix")
        else:
            # Log error and exit program
            logging.error(not_10x_folder_error)
            raise SystemExit(1)
    elif (os.path.isfile(cellranger_multi_folder_bam)) &  (os.path.isdir(cellranger_multi_folder_filtered)):
        cellranger_folder_bam = cellranger_multi_folder_bam
        cellranger_folder_raw = cellranger_multi_folder_filtered
        logging.info("This is a Cell Ranger Multi Directory")
    else:
        logging.error(not_10x_folder_error)
        raise SystemExit(1)

    # Check if the "raw_feature_bc_matrix" directory exists
    if not os.path.isdir(cellranger_folder_raw):
        # Log error and exit 
        logging.error(not_10x_folder_error)
        raise SystemExit(1)
    
    #verify outdir doesn't exist, and if it does create a datetime appended version instead
    if os.path.isdir(args.outdir):
        new_outdir = args.outdir+"_"+datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if os.path.isdir(new_outdir):
            logging.error('Outdir already exists')
            raise SystemExit(1)
        else:
            args.outdir = new_outdir
            logging.info(f"Outdir exists writing to {new_outdir}")
    os.mkdir(args.outdir)

    #miRNA relevant arguments -fixing any that need to be replaced before launching
    miRNA_relevant = {x:i for x,i in args.__dict__.items() if x in ["outdir", "write_bam_output"]}
    miRNA_relevant["matrix_folder"] = cellranger_folder_raw
    miRNA_relevant["outdir"] = miRNA_relevant["outdir"]+"/miRNA"

    #deg relevant 
    deg_relevant = {x:i for x,i in args.__dict__.items() if x in ["outdir", "write_bam_output", "cores", "TSSgtf"]}
    deg_relevant["outdir"] = deg_relevant["outdir"]+"/degradation"
    deg_relevant["matrix_folder"] = cellranger_folder_raw

    if args.miRNAgff3 != None:
        temp_string = "".join([cellranger_folder_bam, " "+ args.miRNAgff3]+[f" --{i} {z}" for i,z in miRNA_relevant.items()]) 
        logging.info(f"Beginning count miRNA with the following arguments: {temp_string}")
        count_miRNAs(temp_string.split(" "))
        #os.system(f"python {pkg_resource_mirna} {cellranger_folder_bam} {args.miRNAgff3} {temp_string}")
    if args.TSSgtf != None:
        temp_string = "".join([cellranger_folder_bam] + [f" --{i} {z}" for i,z in deg_relevant.items()])
        logging.info(f"Beginning deg count with the following arguments: {temp_string}")
        deg_count_with_UMIs(temp_string.split(" ")) 
        #os.system(f"python {pkg_resource_deg} {cellranger_folder_bam} --mtx True {temp_string}")
    elif args.miRNAgff3 == None:
        logging.error("You did not supply a TSS gtf or miRNA gff3 - please supply at least one of these to start Clippings")
        raise SystemExit(1)


if __name__ == "__main__":
     main()
