[![matlab version](	https://img.shields.io/badge/Matlab-2016b-blue.svg)](https://github.com/AswinVasudevan21/CompressiveSensingDictionaryLearning/blob/master/README.md)

# Compressive Sensing Dictionary Learning

### Objective
The dataset we are using originates from a collection of electrophysiological data received from electrode implants in the human brain. The implant in the brain collects data and executes basic processing to allow the device to run properly. The electrode implant sends a wireless signal to a receiver outside the human body. This receiver receives the data from the implant and executes a post-processing procedure on the data collected from the electrode. It is in this external receiver that detection of the neuron spikes and sorting of the dataset takes place.

### Pre Processing
The raw signal received from the implant is extremely noisy and difficult to distinguish. This noise originates from spikes from far away neurons and other random electrical noises. However, compressive sensing can only operate correctly if the signal is sparse. A sparse signal is necessary in order to have a favorable recovery. To create this sparse signal we must apply a fir1 bandpass filter, detect spikes to make the signal sparse, and remove the excess noise by keeping the detected spikes.

### Raw Signal and Cleaned Signal 
<img height="150px" src="https://github.com/AswinVasudevan21/CompressiveSensingDictionaryLearning/blob/master/raw.png">
<img height="150px" src="https://github.com/AswinVasudevan21/CompressiveSensingDictionaryLearning/blob/master/cleaned.png">

### Compressive Sensing:
Compressive sensing theory was developed upon the idea that many signals can be represented using only a few non-zero coefficients in a suitable basis or dictionary. We used L0 and L1 minimization techniques as expalined in our report. 

###  Sampling Methods
The sampling methods we have used for construction of sensing matrix are:
    
    • Random SubSampling : It is achieved by using DCT on the signal and randomly permute
      them to form sensing matrix.
    • Random Gaussian Sampling : Here the randomness is achieved by permutation on the signal
      using Matrix norm.
      
### Conclusion:
We noticed that FISTA algorithm was the fastest and giving consistently good results. We also noticed that Random Subsampling Matrix was giving good results as a sampling method. Gaussian was unpredictable and did not give good results even when values of M were increased to sampling 0.8 of the values.
