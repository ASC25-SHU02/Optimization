[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11046885.svg)](https://doi.org/10.5281/zenodo.11046885)

# KEY POINTS
1. Building indexes is NOT included in elpsed time.
2. Call for software packaging.

# BUILDING
1. download [HISAT-3N](http://daehwankimlab.github.io/hisat2/hisat-3n/), Ver latest
2. download [SRA-Toolkits](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit), Ver3.1.1
3. download [poly dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7051146), change the URL with 46/47/48
    - download with https: `prefetch [SSR]`
    - convert sra to fastq: `fastq-dump [SRR]`
    - convert sra to fasta: `fastq-dump --fasta [SRR]`
4. download [samtools](https://github.com/samtools/samtools/blob/develop/INSTALL) Ver1.21
    - have TINY demand for dependency: `sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev libdeflate-dev`
    - my configuration command: `./configure CC=mpicc`
    - make and make install with `-j $(nproc)`
    > may ADD a `--prefix` to avoid super user privildge.
5. Do accordingly to the handout(tips below)
    - `.fa` in handout is the `.fasta` on this machine
    - `samtools faidx SRR23538290.fasta` get a fai file
6. download [cutseq](https://github.com/y9c/cutseq)
> [!IMPORTANT]
> Others encapsulated in a script... See `./script`


# m<sup>5</sup>C-UBSseq

## Changelog

- 4/23/2024: rewrite code using polars

## workflow

[![](./docs/flow.svg)](https://github.com/y9c/m5C-UBSseq)

## Citation

- cite this software

  ```BibTex
  @software{chang_y_2024_11046885,
      author    = {Chang Y},
      title     = {y9c/m5C-UBSseq: V0.1},
      publisher = {Zenodo},
      version   = {v0.1},
      doi       = {10.5281/zenodo.11046885},
      url       = {https://doi.org/10.5281/zenodo.11046885}
  }
  ```

- cite the method

  ```BibTex
  @article{dai_ultrafast_2024,
      title = {Ultrafast bisulfite sequencing detection of 5-methylcytosine in {DNA} and {RNA}},
      url = {https://www.nature.com/articles/s41587-023-02034-w},
      doi = {10.1038/s41587-023-02034-w},
      author = {Dai, Qing and Ye, Chang and Irkliyenko, Iryna and Wang, Yiding and Sun, Hui-Lung and Gao, Yun and Liu, Yushuai and Beadell, Alana and Perea, José and Goel, Ajay and He, Chuan},
      date = {2024-01-02},
  }
  ```

&nbsp;

<p align="center">
<img
  src="https://raw.githubusercontent.com/y9c/y9c/master/resource/footer_line.svg?sanitize=true"
/>
</p>
<p align="center">
Copyright &copy; 2021-present
<a href="https://github.com/y9c" target="_blank">Chang Y</a>
</p>
<p align="center">
