# PhysiCell-GBM-OV
This project is for the GBM-OV infiltration study. To run, add the code in the repository to the sample_projects folder in 
PhysiCell_1.6.0 and add the following code below to all make files, including "PhysiCell/sample_projects/Makefile-default", "PhysiCell/Makefile-backup" and "PhysiCell/Makefile":

Mitra-biotech-sim:
	cp ./sample_projects/Mitra_biotech_sim/custom_modules/* ./custom_modules/
	touch main.cpp && cp main.cpp main-backup.cpp
	cp ./sample_projects/Mitra_biotech_sim/main-Mitra-biotech-sim.cpp ./main.cpp 
	cp Makefile Makefile-backup
	cp ./sample_projects/Mitra_biotech_sim/Makefile .
	cp ./config/PhysiCell_settings.xml ./config/PhysiCell_settings-backup.xml 
	cp ./sample_projects/Mitra_biotech_sim/config/* ./config/
