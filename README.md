# J Streaming Benchmarks

SC 2020 submission for J Streaming Code.

git-repository: https://github.com/HPCRL/jstream_SC2020_AE

## Building 

J Stream requires Intel C Compiler to compile. Run make command to compile.

```
make
```

## Artifact Evaluation Usage

Scripts require python3. Use diwnload.sh to download input matrix files from suitesparse.

```
./download.sh
```

Use run.sh script to run benchmark for SC results and give the name of folder containing input files as parameter to the script.

```
./run.sh path_to_folder_contatining_mtx_files
```

Before running new tests delete the results from the previous tests by deleting results folder.



