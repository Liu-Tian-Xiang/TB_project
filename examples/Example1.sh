
mv ./BandData ./trash_BandData

#Parameters for the python script to read
<<com 
maxDataTag = 1
minDataTag = 1
com

DataTag=1
mpirun -np 4 ./ecode --run 2222 --dataFlag 00- --dataTag $DataTag --tmpName $DataTag --range -100,400 --BandOffSet 0.3 --inpfile ./Example1_input.txt --AutoMatic 0 --USE_3Dpbc 1 --USE_MaterialMix 0 --USE_FixedSize_or_EnergyRange 0 --USE_OneRange_or_MoreRanges 1 --USE_bulk_states_construction 0 --USE_blocks_transform 0 --USE_blocks_transform_bulk_states 0 --USE_parallel_output 0 --USE_parallel_K_points 0 --USE_distributed 1 --USE_long_parameters 1 --appending

python3 Example1.py

