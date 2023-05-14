# J Streaming SpMM

J Streaming Code for SpMM (sparse-dense matrix-matrix multiplication)

git-repository: https://github.com/HPCRL/jstream_SpMM

## Building 

J Stream requires Intel C Compiler to compile. Run make command to compile.

```
make
```

## Example Usage (for matrices used in SC 2020 paper)

Scripts require python3. Use download.sh to download input matrix files from Suitesparse.

```
./download.sh
```

Use run.sh script to run benchmark for SC results and give the name of folder containing input files as parameter to the script.

```
./run.sh path_to_folder_contatining_mtx_files
```

Before running new tests delete the results from the previous tests by deleting results folder.



