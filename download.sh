cd data
mkdir download
cd download
files="
https://suitesparse-collection-website.herokuapp.com/MM/Um/2cubes_sphere.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/vanHeukelum/cage12.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/Williams/cant.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/Williams/consph.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/Williams/cop20k_A.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/Oberwolfach/filter3D.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/GHS_psdef/hood.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/JGD_Homology/m133-b3.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/Williams/mac_econ_fwd500.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/QLi/majorbasis.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/Williams/mc2depi.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/Um/offshore.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/Pajek/patents_main.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/Williams/pdb1HYS.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/FEMLAB/poisson3Da.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/Boeing/pwtk.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/Bova/rma10.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/Hamm/scircuit.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/DNVS/shipsec1.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/Williams/webbase-1M.tar.gz
https://suitesparse-collection-website.herokuapp.com/MM/SNAP/web-BerkStan.tar.gz
"

for i in $files
do
	wget $i
	
done

for i in ./*.tar.gz;
do
	tar -xzf $i

done

mv */*mtx .


mtx="
2cubes_sphere
cage12
cant
consph
cop20k_A
filter3D
hood
m133-b3
mac_econ_fwd500
majorbasis
mc2depi
offshore
patents_main
pdb1HYS
poisson3Da
pwtk
rma10
scircuit
shipsec1
webbase-1M
web-BerkStan
"

for i in $mtx;
do
	mv $i.mtx ..
done


rm -rf download
cd ..

