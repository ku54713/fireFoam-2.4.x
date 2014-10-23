title = "K25, 97 gpm"

generate:
	generateLookupTableVolumeFlux.py --files k25_97gpm.dat --kFactor 25 --measurementRadius 0.254 --radius 0.254 --aziStep 1 --eleStep 1 --plot --table --out lookup --gpm 97

plot:
	plotLookupTable.py --dataName avgFlux    --legend "Flux (L/s/m2)"     --title $(title) --max 30 --min 0 --image K25-97gpm.jpg

processTube:
	processTubeVolumeFlux.py --file K8-pendent-VolFlux.txt --kFactor 8 --pressure 14 --radius 0.254 --eleStep 2 --aziStep 2 --plot --quadrant 1
plotTube:
	plotLookupTable.py --dataName avgFluxTube   --legend "Tube Flux (L/s/m2)"    --title $(title) --max 11 --min 1 --image K8.jpg
access:
	accessLookupTable.py 
