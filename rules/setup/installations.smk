
rule install_tama:
    output:
        touch(os.path.join(dir.tools_tama,"tama_installed.done"))
    shell:
        """
        if [ -d {dir.tools_tama} ]; then rm -rf {dir.tools_tama}; fi
        git clone https://github.com/GenomeRIK/tama.git {dir.tools_tama} 
        """
