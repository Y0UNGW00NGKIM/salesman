Starter code and data for traveling salesman problem


Files in this directory:

* datareader.cpp : example code to read in the data files (use Makefile)
* datareader.py  : example code to read in the data files
* cities23.dat : list of coordinates for 23 cities in North America
* cities150.dat : 150 cities in North America
* cities1k.dat : 1207 cities in North America
* cities2k.dat : 2063 cities around the world
* routeplot.py : code to plot the globe and salesman's path<br>
usage:<br>
python routeplot.py cities.dat [cities2.dat] -r [="NA"],"World"'<br>
NA = North America, World = Mercator projection of the whole earth
* earth.C : (just for fun) plotting the globe in ROOT


File      Initial     Final       Runtime
cities150 334030.3458 48236.92902 0.137502895
cities1k 2769310.102 132665.0757 0.721354195
cities2k 12003328.79 469594.6222 1.743913952

How to Run:
make clean
make

./sales -i cities150.dat --sweeps 60 --melt 20 --max_seg 200
./sales -i cities1k.dat --sweeps 30 --melt 20 --max_seg 200
./sales -i cities2k.dat --sweeps 25 --melt 20 --max_seg 200

For plots:
python3 routeplot.py cities2k_init.dat cities2k_opt.dat -w --out cities2k.pdf  --initial_km 12003328.79 --final_km 469594.6222 --no_show
python3 plot_schedule.py an2k.dat --out an2k.png --no_show

Where initial_km and final_km are derived from earlier commands. 

Trial configurations are generated using 2 opt move, where I picked two indices along current tour and reverse the intermediate segment limited by --max_seg to keep moves local. Each tour is accepted if it shortens total distance. Otherwise, it's accepted with probability exp(-delta L / T) while temperature is reduced geometrically by T <- alpha * T during annealing