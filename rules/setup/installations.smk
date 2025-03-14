
localrules: install_sqanti, install_tama

rule install_tama:
    output:
        touch(os.path.join(dir.tools_tama,"tama_installed.done"))
    conda:
	f"{dir.envs}/git.yaml"
    shell:
        """
        if [ -d {dir.tools_tama} ]; then rm -rf {dir.tools_tama}; fi
        git clone https://github.com/GenomeRIK/tama.git {dir.tools_tama} 
        """

rule install_sqanti:
    output:
        touch(os.path.join(dir.tools_sqanti,"sqanti_installed.done"))
    conda:
        f"{dir.envs}/git.yaml"    
    shell:
        """
        git clone https://github.com/ConesaLab/SQANTI3.git {dir.tools_sqanti}
        """
