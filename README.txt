This is a README file for the paper with Noam from 2016 to 2019 on analyzing binary halo and their satellite halo population distribution. 

The following is the order and purpose of each python code. Python 3 is used throughout. 

######################## z = 0 CODE ######################## 

----------- script.py

It has a for loop over redshift z, the purpose of which is to identify halo pair of the similar properties in mass and distance to Local Group at each redshift z. Later we abandoned this approach and concentrated on tracing back the identified system from z=0 to higher redshift. We set the processornum to 1 to only process the z=0 data.

## massfilter() function truncate the entire halo catalog into several mass categories. This helps speed up the exhaustive search of the algorithm
## checknosubhalo() get rid of potential subhaloes that are "paired" with either their host or some other haloes
## LGbruteforceFinder() is the MAIN halo pair finder. It checks the box which contains a possible halo LG pair and eliminate those which have either an overmassive halo or a heavy sub halo of LG mass (heavier than the minor partner), or a heavy non-sub halo of LG mass. 
## LGbruteforceFinder() produces file filename_high_res_bf = "LG_HR_bruteforce.dat" which contains the pair indices, which is numbered to be 2998 pairs or 5996 haloes at z=0
Converting the indices to actual halo file is the high_res_trueLG_pop = "high_res_trueLG_pop.dat" file, which contains 5996 lines plus one header line
## after LGbruteforceFinder(), script.py basically has ended. there are some plotting routines for the true LG haloes but are not used. 

----------  nonindexed_smallhaloFinder.py

This is the second part of the LGB halo finder program. The first part we have identified 2998 halo pairs, but some do not satisfy the criteria that there should be no major satellite (>.5 minor host mass). This is to keep the dynamics even clearer by having two distinctly larger haloes. However, there is no way to do this in the last program, because in order to do this we have to first search for large satellites in the system. 

We do this similarly in the last program in LGbruteforceFinder(), first define a box around the halo pair, then check the distances of all the haloes inside the box to the binary center. Then we note the negative data point (system has a large halo within dsep from binary center), and put their line number in the  "high_res_trueLG_pop.dat" file in a separate "neglist.dat". From now on, we uses only "high_res_trueLG_pop.dat" file and "neglist.dat" file and not the hrtrueLG_wosubhalo_biggap_filtered = "high_ res_trueLG_nosubhalo_biggap_filtered_pop.dat" file, because how the halo pair is numbered.

Number of total pairs (including ones with big satellite) (z=0): 2998 pairs or 5996 haloes 
Number of those without big satellites (z=0): 2252 pairs or 4504 haloes 

So for example in the folder of the satellite there should be only 2252 files, each containing the satellite in the system. But the file name which contains the halo pair index is from 0 to 2997, with those with big satellites removed from the folder.

----------  overlapping.py

This program looks for the isolated haloes in the catalog of LG mass range. 

first the program finds everything in the LG mass range halo file (no sub) that is not part of the sample of hrtrueLG_wosubhalo_biggap_filtered. These are the candidates for isolated haloes, they are stored in isolated_lg_file_nosub. 

## truISOLG_finder(): Then we feed them through truISOLG_finder(). truISOLG_finder() first sort through the sample and order the pairs by their separation distance. Then for each sample pair, it goes through the list of potential candidates of isolated haloes, if within 1.5 dsep radius of a candidate there is no halo larger than it, then it satisfy the criteria at the given dsep. If at a given dsep, there is no isolated halo, we make a note of that. This info is stored in isolated_lg_file_at_a_Dsep.

## fakebin_finder(): for a given sample pair, find among the list of candidate haloes that are isolated at the dsep, and find the one which is closest to the halo mass of the sample host. If a candidate has been used already, we add a random rotation to the candidate's satellites. 

## truISOLG_sat_finder() We find the candidates' satellites and rotate them according to what's decided in the last step. 

----------  overlapping_smallhalosAngleDep.py

After identifying LGB and their satellites from the first two steps, and creating 'artificial pair' with no dynamical interaction between them and are isolated. We can now calculate the signal from the satellite distribution and the control signal. 

----------  lopsidedplot.py

This takes hrcosine_distn_file = "high_res_cosine_distn_file.dat" and hrcosine_distn_file_fake = "high_res_cosine_distn_file_fake_1.5dsep.dat" produced in the last step and calculate the histogram of the cosine distribution and other features and dependencies if needed. We notice that somehow the halo finder catalogued also haloes composed of fewer than 20 particles. We therefore made a cut at 20 particles. This reduced the number of satellites (sample size) from  987, 142 to 656, 514. 

plot_cosine_hist_wsigma() densitylinesplot() rtheta_densityplot_all() produce the three plots for the main signal: the cosine signal and comparison between the pos and neg signal, dependence on distance to host, and lastly the minor vs. major host signal

----------- signaloverlapping.py

plots the 2d signal, 2d control and 2d signal/control 

----------- roche.py

Plots the 2D density distribution and signal for various mass ratio of the pair



######################## MERGER TREE CODE ######################## 

----------- merger.py + merger_analysis.py

traces the 2252 pairs back in time. merger_analysis.py also plots the evolution of various features of the pair halo

----------- merger_sat.py + merger_sat_analysis.py

traces the satellites back in time. merger_sat_analysis.py culs the sample so all z have the same size and also cull the sample so all z has mass> 20 particles. 



