#!/usr/bin/env python
import subprocess
import os,sys,yaml
import argparse

def load_configfile(file):
    if not os.path.exists(file):
        sys.exit("Unable to locate the config file; tried %s" % file)
    config = yaml.load(open(file), Loader=yaml.FullLoader)
    return config

def get_snakefile(file="snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

def run_smk(conf,jobs,dryrun,unlock,rerun_incomplete,slurm):
    config = load_configfile(conf)
    sf = get_snakefile()
    conda_prefix = os.path.join(config["required"].get('toolsdir', None),"conda_envs")
    cmd = f'snakemake -s {sf} --configfile {conf} --use-conda --conda-prefix {conda_prefix} -j {jobs}'
    if dryrun:
        cmd += ' -n'
    if unlock:
        cmd += ' --unlock'
    if rerun_incomplete:
        cmd += ' --ri'
    if slurm:
        cmd += ' --slurm --default-resources slurm_account=gge'
    print(f"Running command: {cmd}")
    subprocess.run(cmd, shell=True)
    return

def argparser():
    parser = argparse.ArgumentParser(description='Run snakemake')
    parser.add_argument('--config', '-c', type=str, help='Config file', required=True)
    parser.add_argument('--dryrun', '-n', action='store_true', help='Dry run')
    parser.add_argument('--unlock', '-u', action='store_true', help='Unlock')
    parser.add_argument('--rerun_incomplete', '-R', action='store_true', help='Rerun incomplete')
    parser.add_argument('--slurm', '-s', action='store_true', help='Run on slurm')
    parser.add_argument('--jobs', '-j', type=int, help='Number of jobs',default=8)
    return parser

if __name__ == '__main__':
    parser = argparser()
    args = parser.parse_args()
    run_smk(args.config, args.jobs, args.dryrun, args.unlock, args.rerun_incomplete, args.slurm)