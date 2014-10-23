title = "K14, 120 gpm"

generate:
	generateLookupTableVolumeFlux.py --files k14_115gpm.dat --kFactor 14 --measurementRadius 0.254 --radius 0.254 --aziStep 1 --eleStep 1 --plot --table --out lookup --gpm 120

plot:
	plotLookupTable.py --dataName avgFlux    --legend "Flux (L/s/m2)"     --title $(title) --max 30 --min 0 --image K14.png

processTube:
	processTubeVolumeFlux.py --file K8-pendent-VolFlux.txt --kFactor 8 --pressure 14 --radius 0.254 --eleStep 2 --aziStep 2 --plot --quadrant 1
plotTube:
	plotLookupTable.py --dataName avgFluxTube   --legend "Tube Flux (L/s/m2)"    --title $(title) --max 11 --min 1 --image K8.jpg
access:
	accessLookupTable.py 
