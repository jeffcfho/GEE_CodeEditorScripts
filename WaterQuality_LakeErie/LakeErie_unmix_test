/* Testing the GEE Spectral Unmixing algorithm on an algal bloom image */

var training = ee.Image('LT5_L1T_TOA/LT50190312011246EDC00')
Map.addLayer(training,{bands:['B3','B2','B1'], 'gamma': 2},'True Color')
Map.centerObject(training)

// Use the inspector to get some sample endmembers.
var bloom = [0.11526162177324295, 
             0.10372377187013626,
             0.06162901595234871,
             0.029143797233700752,
             0.002098856261000037,
             295.09185791015625,
             -0.0009386111050844193];
var water = [0.10245747864246368, 
             0.0703064575791359,
             0.044728368520736694,
             0.025622442364692688,
             0.004366065841168165,
             293.7694396972656,
             -0.0009386111050844193];
var sediment = [0.1392693817615509, 
             0.1404828131198883,
             0.11514773219823837,
             0.057314638048410416,
             0.004366065841168165,
             295.529541015625,
             -0.0009386111050844193];

// Now use spectral unmixing to delineate these three types cleanly. 
var unmixed = training.unmix([sediment, bloom, water]);
unmixed = unmixed.rename(['sediment', 'bloom', 'water']);

Map.addLayer(unmixed, {}, 'unmixed');


