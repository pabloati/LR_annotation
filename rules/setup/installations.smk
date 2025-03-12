
rule install_tama:
    output:
        os.path.join(dir.tools_tama,"tama_installed.done")
    shell:
        """
        git clone https://github.com/GenomeRIK/tama.git {dir.tools_tama}
        touch {output}
        """
